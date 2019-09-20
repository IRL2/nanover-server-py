# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.

"""
ASE calculator that implement a wall around the simulation box.

The wall avoids particle to fly out of the box in non-periodic system. Instead,
the particles bounced against the wall.
"""

from collections import Optional
from ase.calculators.calculator import Calculator, all_changes
from ase import Atoms


class VelocityWallCalculator(Calculator):
    def __init__(
            self,
            calculator: Optional[Calculator] = None,
            atoms: Optional[Atoms] = None,
            **kwargs
    ):
        super().__init__(**kwargs)
        self._calculator = calculator
        self.atoms = atoms

    def calculate(self, atoms: Atoms = None, properties=('energy', 'forces'),
                  system_changes=all_changes):

        energy = 0.0
        if atoms is None:
            atoms = self.atoms
        if atoms is None:
            raise ValueError('No ASE atoms supplied to IMD calculation, and no ASE atoms supplied with initialisation.')

        forces = np.zeros((len(atoms), 3))

        if self.calculator is not None:
            self.calculator.calculate(atoms, properties, system_changes)
            if 'energy' in properties:
                energy = self.calculator.results['energy']
            if 'forces' in properties:
                forces = self.calculator.results['forces']

        if 'energy' in properties:
            self.results['energy'] = energy
        if 'forces' in properties:
            self.results['force'] = forces

        self._bounce_atoms(atoms)

    def _bounce_atoms(self, atoms: Atoms):
        # TODO: fail if 
        box = atoms.cell









