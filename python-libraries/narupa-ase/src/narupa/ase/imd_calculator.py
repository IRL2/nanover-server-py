# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.

"""
Provides an implementation of IMD force field in ASE.
"""
from typing import Optional

import numpy as np
from ase import Atoms
from ase.calculators.calculator import Calculator

from narupa.imd.imd_force import calculate_imd_force
from narupa.imd.imd_service import ImdService
import narupa.ase.converter as converter


class ImdCalculator(Calculator):
    """
    Implementation of IMD as an ASE calculator.

    Given another ASE calculator to act as the internal calculator, will compute the external energy
    and forces via the IMD service, and add them to the internal force calculations.

    :param imd_service: The IMD service from which to retrieve interactions to apply as interactive forces.
    :param calculator: An existing ASE calculator to perform internal energy calculation.
    :param atoms: An ASE atoms object to use.
    :param kwargs: Key word args passed to the base calculator.

    """

    def __init__(self, imd_service: ImdService, calculator: Optional[Calculator] = None, atoms: Optional[Atoms] = None,
                 **kwargs):
        """


        """
        super().__init__(**kwargs)
        self._service = imd_service
        self.atoms = atoms
        self._calculator = calculator
        self.implemented_properties = ('energy', 'forces')

    @property
    def calculator(self) -> Calculator:
        """
        The internal ASE calculator being used.
        :return: ASE calculator being used to compute internal forces.
        """
        return self._calculator

    def calculate(self, atoms: Atoms = None, properties=('energy', 'forces'),
                  system_changes=None):

        energy = 0.0
        forces = np.zeros((len(atoms), 3))

        if self.calculator is not None:
            self.calculator.calculate(atoms, properties, system_changes)
            if 'energy' in properties:
                energy = self.calculator.results['energy']
            if 'forces' in properties:
                forces = self.calculator.results['forces']

        imd_energy, imd_forces = self._calculate_imd(atoms)
        if 'energy' in properties:
            self.results['energy'] = energy + imd_energy
        if 'forces' in properties:
            self.results['forces'] = forces + imd_forces

    def _calculate_imd(self, atoms):
        if atoms is None:
            atoms = self.atoms
        if atoms is None:
            raise ValueError('No ASE atoms supplied to IMD calculation, and no ASE atoms supplied with initialisation.')

        # convert positions to the one true unit of distance, nanometers.
        positions = atoms.positions * converter.AngToNm
        # masses are in amu, which is fine.
        masses = atoms.get_masses()

        energy_kjmol, forces_kjmol = calculate_imd_force(positions, masses, self._service.interactions.values())
        ev_per_kjmol = 0.01036
        # convert back to ASE units (eV and Angstroms).
        energy = energy_kjmol * ev_per_kjmol
        forces = forces_kjmol * ev_per_kjmol / converter.NmToAng
        return energy, forces
