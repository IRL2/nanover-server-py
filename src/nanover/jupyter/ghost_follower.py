import numpy as np

from nanover.core import AppServerMinimalImd
from nanover.imd import ParticleInteraction
from nanover.jupyter import ImdAgent
from nanover.jupyter.ghosts import GhostMoleculeObject
from nanover.jupyter.utilities import TransformsUtility
from nanover.trajectory import FrameData


class GhostFollowerAgent(ImdAgent):
    key: str
    ghost: GhostMoleculeObject

    def __init__(self, app_server: AppServerMinimalImd):
        super().__init__(app_server)
        self.transforms = TransformsUtility(app_server)

    def setup(self, *, key: str, ghost: GhostMoleculeObject):
        self.key = key
        self.ghost = ghost

    def update_interactions(self, full_frame: FrameData, frame_update: FrameData):
        key = self.key

        ghost_data = self.ghost.ghost_data
        atom_indices = ghost_data.system_atom_indices

        # target positions are original ghost positions transformed by ghost transform
        target_positions = self.transforms.fetch_transform(
            self.ghost.key
        ).points_local_to_parent(ghost_data.ghost_atom_positions)
        real_positions = full_frame.particle_positions[atom_indices]

        target_centroid = target_positions.mean(axis=0)
        real_centroid = real_positions.mean(axis=0)

        self.ghost.visuals.update_line(
            f"{key}.follow.centroid",
            positions=[real_centroid, target_centroid],
            size=0.01,
            color=[1.0, 0, 0, 1.0],
        )
        self.interactions.update_interaction(
            f"{key}.follow.centroid",
            ParticleInteraction(
                position=target_centroid,
                particles=[int(x) for x in atom_indices],
                type="spring",
                scale=500,
                max_force=100,
            ),
        )

        # rotational following if centroid is close enough
        close = np.linalg.norm(real_centroid - target_centroid, axis=0) < 1

        # find target positions ignoring centroid differences
        rotational_target = target_positions - target_centroid
        rotational_real = real_positions - real_centroid
        rotational = real_positions + (rotational_target - rotational_real)

        for i, index in enumerate(atom_indices):
            if close:
                self.ghost.visuals.update_line(
                    f"{key}.follow.{i}",
                    positions=[rotational[i], real_positions[i]],
                    size=0.01,
                    color=[1.0, 0, 0, 1.0],
                )
                self.interactions.update_interaction(
                    f"{key}.follow.{i}",
                    ParticleInteraction(
                        position=rotational[i],
                        particles=[int(index)],
                        type="spring",
                        scale=100,
                        max_force=50,
                    ),
                )
            else:
                self.ghost.visuals.remove_line(f"{key}.follow.{i}")
                self.interactions.remove_interaction(f"{key}.follow.{i}")
