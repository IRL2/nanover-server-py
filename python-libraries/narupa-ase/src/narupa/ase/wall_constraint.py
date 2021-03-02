import numpy as np
from ase.atoms import Atoms
from ase.cell import Cell


class VelocityWallConstraint:
    """
    Reflective walls implemented as a constraint.
    """

    def __init__(self):
        pass

    def adjust_momenta(self, atoms: Atoms, momenta: np.ndarray):
        box = atoms.cell
        self._validate_box(box)
        positions = atoms.get_positions()
        for dimension in range(3):
            box_max = box[dimension][dimension]
            left = np.logical_and(positions[:, dimension] <= 0,
                                  momenta[:, dimension] < 0)
            right = np.logical_and(positions[:, dimension] >= box_max,
                                   momenta[:, dimension] > 0)
            mask = np.logical_or(left, right)
            momenta[mask, dimension] *= -1

    def adjust_positions(self, atoms: Atoms, positions: np.ndarray):
        pass

    def adjust_forces(self, atoms: Atoms, forces: np.ndarray):
        pass

    @staticmethod
    def _validate_box(cell: Cell):
        """
        Raise an exception if the box is not compatible with the walls.
        """
        if np.isclose(cell.volume, 0):
            raise ValueError('The simulation box has a null volume.')
        if not np.allclose(cell.angles(), [90, 90, 90]):
            raise ValueError('VelocityWall is only compatible with orthorhombic boxes.')
