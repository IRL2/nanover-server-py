"""
Provides an implementation of IMD force field in ASE.
"""
import numpy as np
from ase import Atoms
from ase.calculators.calculator import Calculator

from narupa.imd.imd_force import calculate_imd_force
from narupa.imd.imd_service import ImdService
import narupa.ase.converter as converter


class ImdCalculator(Calculator):
    """
    Implementation of IMD as an ASE calculator.

    """

    def __init__(self, imd_service: ImdService, calculator: Calculator=None, atoms=None, **kwargs):
        super().__init__(**kwargs)
        self._service = imd_service
        self.atoms = atoms
        self._calculator = calculator

    @property
    def calculator(self):
        return self._calculator

    def calculate(self, atoms:Atoms=None, properties=('energy', 'forces'),
                  system_changes=None):

        energy = 0.0
        forces = np.zeros((len(atoms), 3))

        if self.calculator is not None:
            self.calculator.calculate(atoms, properties, system_changes)
            energy = self.calculator.results['energy']
            forces = self.calculator.results['forces']

        imd_energy, imd_forces = self.calculate_imd(atoms)
        self.results['energy'] = energy + imd_energy
        self.results['forces'] = forces + imd_forces

    def calculate_imd(self, atoms):
        if atoms == None:
            atoms = self.atoms
        if atoms == None:
            raise ValueError('No atoms supplied to IMD calculation, and no atoms supplied with initialisation.')

        # convert positions to the one true unit, nanometers.
        positions = atoms.positions * converter.AngToNm
        # masses are in amu, which is fine.
        masses = atoms.get_masses()

        energy_kjmol, forces_kjmol = calculate_imd_force(positions, masses, self._service.interactions.values)
        ev_per_kjmol = 0.01036
        # convert back to ASE units (eV and Angstroms).
        energy = energy_kjmol * ev_per_kjmol
        forces = forces_kjmol * ev_per_kjmol / converter.NmToAng





