from collections.abc import Sequence
from dataclasses import dataclass

import MDAnalysis as mda
import numpy as np
import numpy.typing as npt
from nanover.jupyter import NanoverJupyterUtilities, SceneObjectsUtility
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
        *,
        utilities: NanoverJupyterUtilities,
        key: str,
        atoms: mda.AtomGroup,
        positions: npt.NDArray | None = None,
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
            parent="simulation",
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
