import numpy as np
import numpy.typing as npt

from nanover.imd import ParticleInteraction
from nanover.trajectory import FrameData
from nanover.jupyter import ImdAgent, SceneObjectsUtility
from keyframes import KeyFrame


def fit_keyframe_to_frame(keyframe: KeyFrame, frame: FrameData):
    current = np.array(
        [
            np.mean(frame.particle_positions[target.particles], axis=0)
            for target in keyframe.targets
        ]
    )

    targets = fit_template_points_to_observed(keyframe.centroids, current)

    return current, targets


class SmearAgent(ImdAgent):
    speed = 0.1
    keyframe: KeyFrame | None = None

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.objects = SceneObjectsUtility(self._app_server)

    def close(self):
        super().close()
        self.objects.clear()

    def set_keyframe(self, keyframe: KeyFrame):
        self.interactions.clear()
        self.objects.clear()
        self.keyframe = keyframe

    def update_interactions(self, full_frame: FrameData, frame_update: FrameData):
        if self.keyframe is None:
            return

        # fit keyframe targets to actual positions
        prev_centroids, next_centroids = fit_keyframe_to_frame(
            self.keyframe, full_frame
        )

        # find necessary motions
        deltas = next_centroids - prev_centroids

        # cap motions by speed
        lengths = np.linalg.norm(deltas, axis=1).reshape(-1, 1)
        cappeds = deltas / lengths
        np.clip(lengths, max=self.speed, out=lengths)
        cappeds *= lengths

        # determine final interaction positions
        target_centroids = prev_centroids + cappeds

        with (
            self.interactions as interactions,
            self.objects as objects,
        ):
            points = [target_centroids[i] for i in range(len(self.keyframe.targets))]
            white = [1.0, 1.0, 1.0, 1.0]

            objects.update_line("backbone", positions=points, color=white, size=0.05)

            for i, target in enumerate(self.keyframe.targets):
                error = min(1.0, np.linalg.norm(deltas, axis=1)[i])
                label = f"error: {error:.2g}"
                color = [1.0, 1.0 - error, 1.0 - error, 1.0]
                position = points[i]

                objects.update_shape(i, shape="sphere", position=position, color=color)
                objects.update_label(i, text=label, position=position, color=color)

                if lengths[i] < 0.0001:
                    continue

                interaction = ParticleInteraction(
                    particles=target.particles,
                    position=list(target_centroids[i]),
                    interaction_type="spring",
                    scale=300,
                    max_force=500,
                )
                interactions.update_interaction(
                    f"interaction.REPLAYER.{i}", interaction
                )


def fit_template_points_to_observed(
    template_points: npt.NDArray,
    observed_points: npt.NDArray,
):
    # https://web.stanford.edu/class/cs273/refs/umeyama.pdf
    # https://doi.org/10.1109/34.88573
    p = template_points.T
    q = observed_points.T

    # centroids
    cen_P = np.mean(p, axis=1).reshape(-1, 1)
    cen_Q = np.mean(q, axis=1).reshape(-1, 1)

    # centered vectors
    X = p - cen_P
    Y = q - cen_Q

    # svd
    U, sigma, Vt = np.linalg.svd(X @ Y.T)

    # reflection correction
    d = np.identity(Vt.T.shape[1])
    d[-1, -1] = np.linalg.det(Vt.T @ U.T)

    # rotation and translation
    R = Vt.T @ d @ U.T
    t = cen_Q - R @ cen_P

    return (R @ p + t).T
