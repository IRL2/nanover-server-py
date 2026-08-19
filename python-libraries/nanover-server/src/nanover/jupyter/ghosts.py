from collections.abc import Iterable, Sequence
from dataclasses import dataclass

import MDAnalysis as mda
import numpy as np
import numpy.typing as npt
from nanover.jupyter import NanoverJupyterUtilities, SceneObjectsUtility
from nanover.trajectory import FrameData
from nanover.utilities.transforms import Transform


@dataclass(kw_only=True)
class GhostMoleculeData:
    system_atom_indices: Sequence[int]
    ghost_atom_positions: npt.NDArray
    ghost_bond_pairs: Sequence[Sequence[int]]

    @classmethod
    def from_frame_data(
        cls, frame_data: FrameData, *, atom_indices: Iterable[int] | None = None
    ):
        if atom_indices is None:
            atom_indices = range(frame_data.particle_count)
        atom_indices = list(atom_indices)

        atom_to_index = {atom: index for index, atom in enumerate(atom_indices)}
        ghost_atom_positions = frame_data.particle_positions[atom_indices]
        ghost_bond_pairs = [
            (atom_to_index[a], atom_to_index[b])
            for a, b in frame_data.bond_pairs
            if a in atom_to_index and b in atom_to_index
        ]

        return cls(
            system_atom_indices=atom_indices,
            ghost_atom_positions=ghost_atom_positions,
            ghost_bond_pairs=ghost_bond_pairs,
        )

    @classmethod
    def from_atom_group(cls, atom_group: mda.AtomGroup):
        ghost_universe = mda.Merge(atom_group)
        return cls(
            system_atom_indices=atom_group.atoms.indices,
            ghost_atom_positions=ghost_universe.atoms.positions / 10,  # angstrom -> nm
            ghost_bond_pairs=ghost_universe.bonds.indices,
        )

    def draw(
        self,
        key,
        *,
        visuals: SceneObjectsUtility,
        parent="simulation",
    ):
        draw_ghost(
            key=key,
            visuals=visuals,
            positions=self.ghost_atom_positions,
            bond_pairs=self.ghost_bond_pairs,
            parent=parent,
        )


@dataclass
class GhostMoleculeObject:
    key: str
    ghost_data: GhostMoleculeData
    utilities: NanoverJupyterUtilities
    visuals: SceneObjectsUtility

    @classmethod
    def from_ghost_data(
        cls,
        key: str,
        *,
        utilities,
        ghost_data: GhostMoleculeData,
        parent="simulation",
    ):
        prefix = f"ghost.{key}"

        # normalise around centroid, determine bounding radius
        centroid = np.mean(ghost_data.ghost_atom_positions, axis=0)
        ghost_data.ghost_atom_positions -= centroid
        radius = np.linalg.norm(ghost_data.ghost_atom_positions, axis=0).max()

        # transform + handle for manipulating it
        utilities.transforms.update_transform(
            prefix,
            transform=Transform.from_translation(centroid),
            parent=parent,
        )
        utilities.handles.update_handle(
            prefix, parent=prefix, sphere=((0, 0, 0), radius)
        )

        visuals = SceneObjectsUtility.from_runner(utilities.runner)

        ghost = cls(
            key=prefix,
            ghost_data=ghost_data,
            utilities=utilities,
            visuals=visuals,
        )
        ghost.redraw()
        return ghost

    def redraw(self):
        self.ghost_data.draw(self.key, visuals=self.visuals, parent=self.key)

    def clear(self):
        self.utilities.transforms.remove_transform(self.key)
        self.utilities.handles.remove_handle(self.key)
        self.visuals.clear()


def draw_ghost(
    key: str,
    *,
    visuals: SceneObjectsUtility,
    positions: npt.NDArray,
    bond_pairs: Sequence[Sequence[int]],
    parent="simulation",
):
    for i, position in enumerate(positions):
        visuals.update_shape(
            f"{key}.{i}",
            position=position,
            size=0.1,
            color=[1.0, 1.0, 1.0, 0.5],
            parent=parent,
        )
    for i, (a, b) in enumerate(bond_pairs):
        visuals.update_line(
            f"{key}.{i}",
            positions=positions[[a, b]],
            size=0.05,
            color=[1.0, 1.0, 1.0, 0.5],
            parent=parent,
        )
