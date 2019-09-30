# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.

"""
ASE calculator that implement a wall around the simulation box.

The wall prevents a particle from flying out of the box in non-periodic system. Instead,
the particles bounce against the wall, preserving velocity.
"""

from typing import Optional, Any
import numpy as np
from ase.calculators.calculator import Calculator, all_changes
from ase import Atoms
from ase.cell import Cell


class VelocityWallCalculator(Calculator):
    def __init__(
            self,
            calculator: Optional[Calculator] = None,
            atoms: Optional[Atoms] = None,
            **kwargs
    ):
        self._parent_calculator = calculator
        super().__init__(**kwargs)
        self.implemented_properties = ('energy', 'forces')
        if self._parent_calculator is not None:
            self.implemented_properties = tuple(
                set(self.implemented_properties)
                | set(self._parent_calculator.implemented_properties)
            )

    def calculate(self, atoms: Atoms, properties=('energy', 'forces'),
                  system_changes=all_changes):
        if atoms is None:
            raise ValueError(
                'No ASE atoms supplied to IMD calculation, '
                'and no ASE atoms supplied with initialisation.'
            )

        energy = 0.0
        forces = np.zeros((len(atoms), 3))
        if self._parent_calculator is not None:
            self._parent_calculator.calculate(atoms, properties, system_changes)
            for key, value in self._parent_calculator.results.items():
                self.results[key] = value
        self.results['energy'] = self.results.get('energy', energy)
        self.results['forces'] = self.results.get('forces', forces)

        self._validate_box(atoms.cell)
        self._bounce_atoms(atoms)

    @staticmethod
    def _bounce_atoms(atoms: Atoms):
        box = atoms.cell
        positions = atoms.get_positions()
        velocities = atoms.get_velocities()
        for dimension in range(3):
            box_max = box[dimension][dimension]
            left = np.logical_and(positions[:, dimension] <= 0,
                                  velocities[:, dimension] < 0)
            right = np.logical_and(positions[:, dimension] >= box_max,
                                   velocities[:, dimension] > 0)
            mask = np.logical_or(left, right)
            velocities[mask, dimension] *= -1
        atoms.set_velocities(velocities)

    @staticmethod
    def _validate_box(cell: Cell):
        """
        Raise an exception if the box is not compatible with the walls.
        """
        if np.isclose(cell.volume, 0):
            raise ValueError('The simulation box has a null volume.')
        if not np.allclose(cell.angles(), [90, 90, 90]):
            raise ValueError('VelocityWall only works for orthorhombic boxes.')
    
    def __getattr__(self, name: str) -> Any:
        """
        Give direct access to the attribute of the parent calculator if they
        are not redefined in this one.

        This is required as some uses expect some specific methods or attributes
        from a specific calculator, and they must be available even if the
        calculator is chained (e.g. OpenMMCalculator.topology).
        """
        return getattr(self._parent_calculator, name)
