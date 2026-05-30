import traceback
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field

import numpy as np
import numpy.typing as npt

from nanover.app import OmniRunner
from nanover.imd import ParticleInteraction
from nanover.utilities.cli import CancellationToken


@dataclass(kw_only=True, eq=False)
class TargetGroup:
    particles: list[int] = field(default_factory=list)
    centroid: npt.NDArray[np.float32]


@dataclass(kw_only=True, eq=False)
class KeyFrame:
    targets: list[TargetGroup] = field(default_factory=list)
    centroids: npt.NDArray[np.float32] = field(init=False)

    def __post_init__(self):
        self.centroids = np.array([target.centroid for target in self.targets])


class SmearAgent:
    @classmethod
    def from_runner(cls, runner: OmniRunner):
        return cls(runner)

    def __init__(self, runner: OmniRunner):
        self._runner = runner
        self._threads = ThreadPoolExecutor(max_workers=1)
        self._cancellation = CancellationToken()
        self._task = None
        self._interactions: set[str] = set()

    def start(self, keyframe: KeyFrame, speed: float, output):
        if self._task is not None:
            return

        publisher = self._runner.app_server.frame_publisher
        imd = self._runner.app_server.imd
        stream = publisher.subscribe_latest_frames(
            frame_interval=0, cancellation=self._cancellation
        )

        def run():
            try:
                for frame in stream:
                    prev_centroids = np.array([
                        np.average(frame.particle_positions[target.particles], axis=0)
                        for target in keyframe.targets
                    ])

                    # fit keyframe targets to actual positions
                    next_centroids = fit_template_points_to_observed(keyframe.centroids, prev_centroids)

                    # find necessary motions
                    deltas = next_centroids - prev_centroids
                    error = np.sum(deltas)

                    # cap motions by speed
                    lengths = np.array([np.linalg.norm(deltas, axis=1)]).transpose()
                    cappeds = deltas / lengths
                    np.clip(lengths, max=speed, out=lengths)
                    cappeds *= lengths

                    # determine final interaction positions
                    target_centroids = prev_centroids + cappeds

                    for i, target in enumerate(keyframe.targets):
                        key = f"interaction.REPLAYER.{i}"
                        self._interactions.add(key)

                        if lengths[i] > 0.0001:
                            imd.insert_interaction(
                                key,
                                ParticleInteraction(
                                    particles=target.particles,
                                    position=list(target_centroids[i]),
                                    interaction_type="spring",
                                    scale=500,
                                    max_force=500,
                                ),
                            )
            except Exception as e:
                with output:
                    print(traceback.print_exc())

        self._task = self._threads.submit(run)

    def stop(self):
        self._cancellation.cancel()
        self._threads.shutdown(wait=True)
        for key in self._interactions:
            self._runner.app_server.imd.remove_interaction(key)


def fit_template_points_to_observed(template_points: npt.NDArray, observed_points: npt.NDArray):
    # minimise (Rp+t - q)**2 -- fit template points onto observed points
    p = template_points.transpose()
    q = observed_points.transpose()

    # centroids
    cen_P = np.mean(p, axis=1).reshape(-1, 1)
    cen_Q = np.mean(q, axis=1).reshape(-1, 1)

    # positions relative to centroids ("centered vectors")
    X = p - cen_P
    Y = q - cen_Q

    # covariance matrix
    S = X @ Y.T
    U, sigma, Vt = np.linalg.svd(S)

    # diagonal matrix
    d = np.eye(Vt.T.shape[1])
    d[-1, -1] = np.linalg.det(Vt.T @ U.T)

    # rotation and translation
    R = Vt.T @ d @ U.T
    t = cen_Q - R @ cen_P

    return (R @ p + t).transpose()
