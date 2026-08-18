from collections.abc import Iterable, Sequence
from dataclasses import dataclass

import MDAnalysis as mda
import numpy as np
import numpy.typing as npt
from nanover.jupyter import NanoverJupyterUtilities, SceneObjectsUtility
from nanover.trajectory import FrameData
from nanover.utilities.transforms import Transform


@dataclass
class GhostMolecule:
    key: str
    positions: npt.NDArray
    atom_indices: Sequence[int]  # global indices
    bond_pairs: Sequence[Sequence[int]]  # ghost indices
    utilities: NanoverJupyterUtilities
    visuals: SceneObjectsUtility

    @classmethod
    def from_mdanalysis(
        cls,
        key: str,
        *,
        utilities: NanoverJupyterUtilities,
        atoms: mda.AtomGroup,
        positions: npt.NDArray | None = None,
        parent="simulation",
    ):
        prefix = f"ghost.{key}"

        # extract ghost molecule
        ghost_universe = mda.Merge(atoms)
        ghost_positions = (
            positions if positions is not None else ghost_universe.atoms.positions / 10
        )  # angstrom -> nm

        # normalise around centroid, determine bounding radius
        centroid = np.mean(ghost_positions, axis=0)
        ghost_positions -= centroid
        radius = np.linalg.norm(ghost_positions, axis=0).max()

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
            positions=ghost_positions,
            atom_indices=atoms.indices,
            bond_pairs=ghost_universe.bonds.indices,
            utilities=utilities,
            visuals=visuals,
        )
        ghost.redraw()
        return ghost

    def redraw(self):
        draw_ghost(
            self.key,
            visuals=self.visuals,
            positions=self.positions,
            bond_pairs=self.bond_pairs,
            parent=self.key,
        )

    def clear(self):
        self.utilities.transforms.remove_transform(self.key)
        self.utilities.handles.remove_handle(self.key)
        self.visuals.clear()


def draw_ghost_from_atom_group(
    key: str,
    *,
    visuals: SceneObjectsUtility,
    atom_group: mda.AtomGroup,
    parent="simulation",
):
    universe = mda.Merge(atom_group)
    positions = universe.atoms.positions / 10
    bond_pairs = universe.bonds.indices

    draw_ghost(
        key,
        visuals=visuals,
        positions=positions,
        bond_pairs=bond_pairs,
        parent=parent,
    )


def draw_ghost_from_frame_data(
    key: str,
    *,
    visuals: SceneObjectsUtility,
    frame_data: FrameData,
    atoms: Iterable[int],
    parent="simulation",
):
    atom_to_index = {atom: index for index, atom in enumerate(atoms)}
    bond_pairs = [
        (atom_to_index[a], atom_to_index[b])
        for a, b in frame_data.bond_pairs
        if a in atom_to_index and b in atom_to_index
    ]

    draw_ghost(
        key,
        visuals=visuals,
        positions=frame_data.particle_positions[np.asarray(atoms)],
        bond_pairs=bond_pairs,
        parent=parent,
    )


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
