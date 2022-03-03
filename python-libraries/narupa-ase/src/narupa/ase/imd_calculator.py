# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.

"""
Provides an implementation of IMD force field in ASE.
"""
import math
from typing import Optional, Dict, Set, Collection

import numpy as np
from ase import Atoms, units  # type: ignore
from ase.calculators.calculator import Calculator, all_changes
from ase.md.md import MolecularDynamics
from ase.md.velocitydistribution import _maxwellboltzmanndistribution
from narupa.imd.imd_force import calculate_imd_force
from narupa.imd.imd_state import ImdStateWrapper
from narupa.imd.particle_interaction import ParticleInteraction
from narupa.trajectory.frame_data import MissingDataError

from . import converter


class ImdCalculator(Calculator):
    """
    Implementation of IMD as an ASE calculator.

    Given another ASE calculator to act as the internal calculator, will compute the external energy
    and forces via the IMD service, and add them to the internal force calculations.

    :param imd_state: A wrapper that provides access to the active interactive forces.
    :param calculator: An existing ASE calculator to perform internal energy calculation.
    :param atoms: An ASE atoms object to use.
    :param dynamics: An ASE dynamics object from which to draw the equilibrium temperature for resetting velocities
    :param reset_scale: A scale factor to apply to velocities when reset.
    :param kwargs: Key word args passed to the base calculator.

    .. seealso::

        The :class:`ImdServer` class makes use of this class, and makes
        running an interactive molecular dynamics simulation in ASE straightforward.
    
    """

    def __init__(self, imd_state: ImdStateWrapper,
                 calculator: Optional[Calculator] = None,
                 atoms: Optional[Atoms] = None,
                 dynamics: Optional[MolecularDynamics] = None,
                 reset_scale=0.5,
                 **kwargs):

        super().__init__(**kwargs)
        self._imd_state = imd_state
        self.atoms = atoms
        self._calculator = calculator
        self.implemented_properties = ('energy', 'forces', 'interactive_energy', 'interactive_forces')
        self._dynamics = dynamics
        self.reset_scale = reset_scale
        self._custom_temperature = None
        self._initialise_velocity_reset()

    @property
    def temperature(self) -> float:
        """
        Gets the temperature used for reinitialising velocities after an interaction.

        By default, it will attempt to use the temperature of the dynamics.
        If a custom temperature has been set by this attributes setter, then that will be used.

        :return: The temperature used for reinitialising velocities after an interaction.
        :raises: AttributeError: If no temperature is defined for this calculator, in the case
            that no dynamics object has been passed, or the dynamics object does not implement the
            temperature or 'temp' attribute.
        """
        if self._custom_temperature is not None:
            return self._custom_temperature

        if self._dynamics is None:
            raise MissingDataError(
                "No temperature has been set, and no molecular dynamics object has been passed to the "
                "IMD calculator.")

        try:
            return self._dynamics.temperature
        except AttributeError:
            try:
                return self._dynamics.temp
            except AttributeError:
                raise MissingDataError("No temperature has been set, and the molecular dynamics object does not "
                                       "appear to set a temperature.")

    @temperature.setter
    def temperature(self, value):
        """
        Sets the temperature used for reinitialising velocities after an interaction. Note that
        if this is set, it will be used instead of the temperature that the dynamics is running at.

        :param value: The custom temperature to use.
        """
        self._custom_temperature = value

    @property
    def reset_temperature(self):
        """
        The temperature this calculator will reset the velocities of atoms interacted with to if the interaction
        is set to reset velocities.

        :return: The reset temperature.
        :raises: Attribute error, if not temperature has been defined.
        """
        return self.temperature * self.reset_scale

    @property
    def calculator(self) -> Calculator:
        """
        The internal ASE calculator being used.

        :return: ASE calculator being used to compute internal forces.
        """
        return self._calculator

    @property
    def interactions(self) -> Dict[str, ParticleInteraction]:
        """
        Fetches a copy of the current interactions.
        """

        return self._imd_state.active_interactions

    def calculate(self, atoms: Atoms = None, properties=('energy', 'forces'),
                  system_changes=all_changes):
        """
        Calculates the given properties of the ASE atoms. The internal molecular calculator is called first,
        and then any interactive forces currently being applied to the system are added.

        Results are stored in the results dictionary, as normal.

        :param atoms: Optional :class:`Atoms` object to perform the calculation on. If no atoms is passed,
            the atoms object passed at initialisation are used.
        :param properties: The properties to calculate. The ImdCalculator support 'energy' and 'forces',
            but will pass any other requested properties to the internal atomic calculator.
            See :func:`~Calculator.calculate` for details.
        :param system_changes: List of what has changed since last calculation. See :func:`~Calculator.calculate` for
            details.

        :raises ValueError: If no ASE atoms are supplied to the calculation, and no ASE atoms were supplied during
            initialisation.
        """
        energy = 0.0
        if atoms is None:
            atoms = self.atoms
        if atoms is None:
            raise ValueError('No ASE atoms supplied to IMD calculation, and no ASE atoms supplied with initialisation.')

        forces = np.zeros((len(atoms), 3))

        if self.calculator is not None:
            self.calculator.calculate(atoms, properties, system_changes)
            energy = self.calculator.results['energy']
            forces = self.calculator.results['forces']

        imd_energy, imd_forces = self._calculate_imd(atoms)
        self.results['energy'] = energy + imd_energy
        self.results['forces'] = forces + imd_forces
        self.results['interactive_energy'] = imd_energy
        self.results['interactive_forces'] = imd_forces

    def _calculate_imd(self, atoms):

        interactions = self.interactions

        self._reset_velocities(atoms, interactions, self._previous_interactions)

        # convert positions to the one true unit of distance, nanometers.
        positions = atoms.positions * converter.ANG_TO_NM
        # masses are in amu, which is fine.
        masses = atoms.get_masses()

        periodic_box_lengths = get_periodic_box_lengths(atoms)
        energy_kjmol, forces_kjmol = calculate_imd_force(positions, masses, interactions.values(),
                                                         periodic_box_lengths=periodic_box_lengths)
        ev_per_kjmol = converter.KJMOL_TO_EV
        # convert back to ASE units (eV and Angstroms).
        energy = energy_kjmol * ev_per_kjmol
        forces = forces_kjmol * ev_per_kjmol / converter.NM_TO_ANG

        # update previous interactions for next step.
        self._previous_interactions = dict(interactions)

        return energy, forces

    def _reset_velocities(self, atoms, interactions, previous_interactions):

        cancelled_interactions = _get_cancelled_interactions(interactions, previous_interactions)
        atoms_to_reset = _get_atoms_to_reset(cancelled_interactions)
        if len(atoms_to_reset) == 0:
            return
        # If no temperature has been defined, we cannot reinitialise velocities.
        # check for temperature before doing anything, so state doesn't change.
        reset_temperature = self.reset_temperature
        _apply_velocities_reset(atoms, atoms_to_reset, reset_temperature)

    def _initialise_velocity_reset(self):
        try:
            temp = self.temperature
        except MissingDataError:
            self._imd_state.velocity_reset_available = False
        self._imd_state.velocity_reset_available = True
        self._previous_interactions = {}


def get_periodic_box_lengths(atoms: Atoms) -> Optional[np.ndarray]:
    """
    Gets the periodic box lengths of an orthorhombic box, in nm, from an ASE atoms collection, if it exists.

    :param atoms: ASE atoms.
    :return: Array of periodic box lengths if periodic boundaries are in use, ``None`` otherwise.
    """
    if not np.all(atoms.get_pbc()):
        if np.any(atoms.get_pbc()):
            raise NotImplementedError(f'Atoms object has periodic unit cell on only some dimensions, which is not '
                                      f'supported.')
        return None
    lengths_angles = atoms.cell.cellpar()
    lengths = np.array(lengths_angles[0:3])
    angles = np.array(lengths_angles[3:])
    if not np.allclose(angles, (90, 90, 90)):
        raise NotImplementedError(
            f'Atoms object has periodic unit cell that is not orthorhombic, which is not supported!')
    return lengths


def _get_cancelled_interactions(interactions, previous_interactions) -> Dict[object, ParticleInteraction]:
    old_keys = set(previous_interactions.keys())
    cancelled_interactions = old_keys.difference(interactions.keys())
    return {key: previous_interactions[key] for key in cancelled_interactions}


def _get_atoms_to_reset(cancelled_interactions) -> Set[int]:
    atoms_to_reset: Set[int] = set()
    for key, interaction in cancelled_interactions.items():
        if interaction.reset_velocities:
            atoms_to_reset = atoms_to_reset.union(interaction.particles)
    return atoms_to_reset


def _apply_velocities_reset(atoms, atoms_to_reset, temperature):
    atoms_to_reset = np.array(list(atoms_to_reset))

    _reset_selection_to_boltzmann(atoms, atoms_to_reset, temperature)
    # now scale the velocities so the exact target temperature is achieved.
    _scale_momentum_of_selection(atoms, atoms_to_reset, temperature)


def _reset_selection_to_boltzmann(atoms: Atoms, selection: Collection[int], temperature: float):
    # TODO importing a private function here... reimplement?
    reset = _maxwellboltzmanndistribution(atoms[selection].get_masses(), temperature * units.kB)
    _apply_momentum_to_selection(atoms, selection, reset)


def _scale_momentum_of_selection(atoms: Atoms, selection: Collection[int], temperature: float):
    scaled_selection = _get_scaled_momentum(atoms[selection], temperature)
    _apply_momentum_to_selection(atoms, selection, scaled_selection)


def _get_scaled_momentum(atoms: Atoms, temperature):
    current_temp = atoms.get_temperature()
    scale = temperature / current_temp
    return atoms.get_momenta() * math.sqrt(scale)


def _apply_momentum_to_selection(atoms: Atoms, selection, momentum):
    m = atoms.get_momenta()
    m[selection] = momentum
    atoms.set_momenta(m)
