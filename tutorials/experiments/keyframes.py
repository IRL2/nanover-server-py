from dataclasses import dataclass, field

import numpy as np
import numpy.typing as npt
from MDAnalysis import Universe

from nanover.imd.imd_state import dict_to_interaction
from nanover.recording import NanoverRecordingReader
from nanover.utilities.change_buffers import DictionaryChange

RESTRAINT_PREFIX = "interaction.MOVEABLE-RESTRAINT"
CHECKPOINT_KEY = "mark.checkpoint"


@dataclass(kw_only=True, eq=False)
class TargetGroup:
    particles: list[int] = field(default_factory=list)
    position: npt.NDArray[np.float32]
    centroid: npt.NDArray[np.float32]


@dataclass(kw_only=True, eq=False)
class KeyFrame:
    targets: list[TargetGroup] = field(default_factory=list)
    centroids: npt.NDArray[np.float32] = field(init=False)

    def __post_init__(self):
        self.centroids = np.array([target.centroid for target in self.targets])


def extract_keyframes(reader: NanoverRecordingReader):
    keyframes = []

    for event in reader.iter_max():
        if event.next_state_event is None:
            continue

        change = DictionaryChange.from_dict(event.next_state_event.message)

        if CHECKPOINT_KEY not in change.updates:
            continue

        targets = []
        positions = event.next_frame.particle_positions

        for key in event.next_state:
            if key.startswith(RESTRAINT_PREFIX):
                restraint = dict_to_interaction(event.next_state[key])
                particles = [int(i) for i in restraint.particles]
                centroid = np.average(positions[particles], axis=0)

                target = TargetGroup(
                    particles=particles,
                    centroid=centroid,
                    position=np.array(restraint.position),
                )
                targets.append(target)

        if targets:
            keyframes.append(KeyFrame(targets=targets))

    return keyframes


def extract_keyframes_brute(reader: NanoverRecordingReader, universe: Universe):
    keyframes = []

    for event in reader.iter_max():
        if event.next_state_event is None:
            continue

        change = DictionaryChange.from_dict(event.next_state_event.message)

        if CHECKPOINT_KEY not in change.updates:
            continue

        targets = []
        positions = event.next_frame.particle_positions

        for residue in universe.residues:
            particles = [int(i) for i in residue.atoms.indices]
            centroid = np.average(positions[particles], axis=0)

            target = TargetGroup(
                particles=particles,
                centroid=centroid,
                position=centroid,
            )
            targets.append(target)

        if targets:
            keyframes.append(KeyFrame(targets=targets))

    return keyframes
