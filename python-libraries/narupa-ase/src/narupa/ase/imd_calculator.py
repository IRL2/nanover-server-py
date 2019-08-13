# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.

"""
Provides an implementation of IMD force field in ASE.
"""
from typing import Optional, Dict, Tuple

from . import converter
import numpy as np
from ase import Atoms
from ase.calculators.calculator import Calculator, all_changes
from narupa.imd.imd_force import calculate_imd_force
from narupa.imd.imd_service import ImdService
from narupa.imd.particle_interaction import ParticleInteraction


def get_periodic_box_lengths(atoms: Atoms) -> Optional[np.ndarray]:
    """
    Gets the periodic box lengths of an orthorhombic box, in nm, from an ASE atoms collection, if it exists.
    :param atoms: ASE atoms.
    :return: Array of periodic box lengths if periodic boundaries are in use, None otherwise.
    """
    if not np.all(atoms.get_pbc()):
        if np.any(atoms.get_pbc()):
            raise NotImplementedError(f'Atoms object has periodic unit cell on only some dimensions, which is not '
                                      f'supported.')
        return None
    lengths_angles = atoms.get_cell_lengths_and_angles()
    lengths = np.array(lengths_angles[0:3])
    angles = np.array(lengths_angles[3:])
    if not np.allclose(angles, (90, 90, 90)):
        raise NotImplementedError(
            f'Atoms object has periodic unit cell that is not orthorhombic, which is not supported!')
    return lengths


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
        self.implemented_properties = ('energy', 'forces', 'interactive_energy', 'interactive_forces')

    @property
    def calculator(self) -> Calculator:
        """
        The internal ASE calculator being used.

        :return: ASE calculator being used to compute internal forces.
        """
        return self._calculator

    @property
    def interactions(self) -> Dict[Tuple[str, str], ParticleInteraction]:
        """
        Returns a shallow copy of the current interactions.
        """

        return self._service.active_interactions

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

        imd_energy, imd_forces = self._calculate_imd(atoms)
        if 'energy' in properties:
            self.results['energy'] = energy + imd_energy
        if 'forces' in properties:
            self.results['forces'] = forces + imd_forces
        if 'interactive_energy' in properties:
            self.results['interactive_energy'] = imd_energy
        if 'interactive_forces' in properties:
            self.results['interactive_forces'] = imd_forces

    def _calculate_imd(self, atoms):

        # convert positions to the one true unit of distance, nanometers.
        positions = atoms.positions * converter.ANG_TO_NM
        # masses are in amu, which is fine.
        masses = atoms.get_masses()

        periodic_box_lengths = get_periodic_box_lengths(atoms)
        interactions = self.interactions
        energy_kjmol, forces_kjmol = calculate_imd_force(positions, masses, interactions,
                                                         periodic_box_lengths=periodic_box_lengths)
        ev_per_kjmol = converter.KJMOL_TO_EV
        # convert back to ASE units (eV and Angstroms).
        energy = energy_kjmol * ev_per_kjmol
        forces = forces_kjmol * ev_per_kjmol / converter.NM_TO_ANG
        return energy, forces
