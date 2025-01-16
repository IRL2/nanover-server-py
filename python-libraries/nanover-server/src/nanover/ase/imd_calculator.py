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

from nanover.imd.imd_force import calculate_imd_force, get_sparse_forces
from nanover.imd.imd_state import ImdStateWrapper
from nanover.imd.particle_interaction import ParticleInteraction
from nanover.trajectory.frame_data import MissingDataError, FrameData

from . import converter


class ImdForceManager:
    """
    A class that calculates and stores the iMD forces and energies for
    an ASE simulation. This class manages data associated with the iMD interactions
    for the :class:`ImdCalculator`.

    :param imd_state: A wrapper that provides access to the active interactive forces.
    :param atoms: An ASE atoms object to use.
    """

    def __init__(self, imd_state: ImdStateWrapper, atoms: Atoms):
        self.atoms = atoms
        self.imd_state = imd_state

        self.total_user_energy = 0.0
        self.user_forces: np.ndarray = np.zeros(self.atoms.positions.shape)

        self._current_interactions: Dict[str, ParticleInteraction] = {}

    def update_interactions(self):
        """
        Update the iMD interaction energies and forces (in ASE units).
        """
        self._update_forces(self.atoms)

    def add_to_frame_data(self, frame_data: FrameData):
        """
        Add the iMD forces and energy to the frame data.

        :param frame_data: The FrameData object to which the iMD results are appended.
        """
        frame_data.user_energy = self.total_user_energy * converter.EV_TO_KJMOL
        user_sparse_indices, user_sparse_forces = get_sparse_forces(self.user_forces)
        frame_data.user_forces_sparse = user_sparse_forces * (
            converter.EV_TO_KJMOL / converter.ANG_TO_NM
        )
        frame_data.user_forces_index = user_sparse_indices

    def _update_forces(self, atoms):
        """
        Get the forces to apply from the iMD service and communicate them

        :param atoms: an ASE atoms object.
        """
        # Get the current interactions from the iMD service (if any)
        interactions = self.imd_state.active_interactions

        # convert positions to the one true unit of distance, nanometers.
        positions = atoms.positions * converter.ANG_TO_NM
        energy = 0.0
        forces = np.zeros(positions.shape)
        # If there are iMD interactions, calculate their forces and energies
        if interactions:
            energy, forces = self._calculate_imd(atoms, positions, interactions)

        # Store the energy and forces
        self.total_user_energy = energy
        self.user_forces = forces

        # update interactions for next frame interval
        self._current_interactions = dict(interactions)

    def _calculate_imd(
        self,
        atoms,
        positions: np.ndarray,
        interactions: Dict[str, ParticleInteraction],
    ):
        """
        A calculate the iMD forces and energies and convert
        the results into ASE units.

        :param atoms: an ASE atoms object.
        :param positions: the positions of the ASE atoms in NanoVer units.
        :param interactions: a dictionary of the interactions being applied
        to the atoms in the system.
        """
        # masses are in amu, which is fine.
        masses = atoms.get_masses()

        periodic_box_lengths = get_periodic_box_lengths(atoms)
        energy_kjmol, forces_kjmol = calculate_imd_force(
            positions,
            masses,
            interactions.values(),
            periodic_box_lengths=periodic_box_lengths,
        )
        ev_per_kjmol = converter.KJMOL_TO_EV
        # convert back to ASE units (eV and Angstroms).
        energy = energy_kjmol * ev_per_kjmol
        forces = forces_kjmol * ev_per_kjmol / converter.NM_TO_ANG

        return energy, forces


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

    def __init__(
        self,
        imd_state: ImdStateWrapper,
        imd_force_manager: Optional[ImdForceManager] = None,
        calculator: Optional[Calculator] = None,
        atoms: Optional[Atoms] = None,
        dynamics: Optional[MolecularDynamics] = None,
        reset_scale=0.5,
        **kwargs,
    ):
        super().__init__(**kwargs)
        self._imd_state = imd_state
        self._imd_force_manager = imd_force_manager
        self.atoms = atoms
        self._calculator = calculator
        self.implemented_properties = [
            "energy",
            "forces",
            "interactive_energy",
            "interactive_forces",
        ]
        self._dynamics = dynamics
        self.reset_scale = reset_scale
        self._custom_temperature = None
        self._initialise_velocity_reset()
        self._pbc_implemented = False

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
                "IMD calculator."
            )

        # Some, but not all, dynamics define a temperature or temp attribute.
        try:
            return self._dynamics.temperature  # type: ignore
        except AttributeError:
            try:
                return self._dynamics.temp  # type: ignore
            except AttributeError:
                raise MissingDataError(
                    "No temperature has been set, and the molecular dynamics object does not "
                    "appear to set a temperature."
                )

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
    def calculator(self) -> Optional[Calculator]:
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

    def calculate(
        self,
        atoms: Optional[Atoms] = None,
        properties=("energy", "forces"),
        system_changes=all_changes,
    ):
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
            raise ValueError(
                "No ASE atoms supplied to IMD calculation, and no ASE atoms supplied with initialisation."
            )

        # Check whether the periodic boundary conditions defined in the atoms object are implemented
        # for the iMD calculator. If they are (or no pbcs are defined), set _pbc_implemented property
        # to true to avoid repeating this check.
        if not self._pbc_implemented:
            get_periodic_box_lengths(atoms)
            self._pbc_implemented = True

        forces = np.zeros((len(atoms), 3))

        if self.calculator is not None:
            self.calculator.calculate(atoms, properties, system_changes)
            energy = self.calculator.results["energy"]
            forces = self.calculator.results["forces"]

        if self._imd_force_manager is not None:
            # Retrieve iMD energy and forces and add to results
            imd_energy = self._imd_force_manager.total_user_energy
            imd_forces = self._imd_force_manager.user_forces
            self.results["energy"] = energy + imd_energy
            self.results["forces"] = forces + imd_forces
            self.results["interactive_energy"] = imd_energy
            self.results["interactive_forces"] = imd_forces

    def update_interactions(self):
        """
        Update the iMD interaction energies and forces (in ASE units)
        via the ImdForceManager, and subsequently resets the velocities
        (if applicable).
        """
        assert self._imd_force_manager is not None
        prev_interactions = self._imd_force_manager._current_interactions
        self._imd_force_manager.update_interactions()
        next_interactions = self._imd_force_manager._current_interactions
        self._reset_velocities(self.atoms, next_interactions, prev_interactions)

    def add_to_frame_data(self, frame_data: FrameData):
        """
        Add the iMD forces and energy to the frame data via the
        ImdForceManager.

        :param frame_data: The FrameData object to which the iMD results are appended.
        """
        assert self._imd_force_manager is not None
        self._imd_force_manager.add_to_frame_data(frame_data)

    def _reset_velocities(self, atoms, interactions, previous_interactions):
        cancelled_interactions = _get_cancelled_interactions(
            interactions, previous_interactions
        )
        atoms_to_reset = _get_atoms_to_reset(cancelled_interactions)
        if len(atoms_to_reset) == 0:
            return
        # If no temperature has been defined, we cannot reinitialise velocities.
        # check for temperature before doing anything, so state doesn't change.
        reset_temperature = self.reset_temperature
        _apply_velocities_reset(atoms, atoms_to_reset, reset_temperature)

    def _initialise_velocity_reset(self):
        try:
            pass
        except MissingDataError:
            self._imd_state.velocity_reset_available = False
        self._imd_state.velocity_reset_available = True


def get_periodic_box_lengths(atoms: Atoms) -> Optional[np.ndarray]:
    """
    Gets the periodic box lengths of an orthorhombic box, in nm, from an ASE atoms collection, if it exists.

    :param atoms: ASE atoms.
    :return: Array of periodic box lengths if periodic boundaries are in use, ``None`` otherwise.
    """
    if not np.all(atoms.get_pbc()):
        if np.any(atoms.get_pbc()):
            raise NotImplementedError(
                "Atoms object has periodic unit cell on only some dimensions, which is not "
                "supported."
            )
        return None
    lengths_angles = atoms.cell.cellpar()
    lengths = np.array(lengths_angles[0:3])
    angles = np.array(lengths_angles[3:])
    if not np.allclose(angles, (90, 90, 90)):
        raise NotImplementedError(
            "Atoms object has periodic unit cell that is not orthorhombic, which is not supported!"
        )
    return lengths


def _get_cancelled_interactions(
    interactions, previous_interactions
) -> Dict[object, ParticleInteraction]:
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


def _reset_selection_to_boltzmann(
    atoms: Atoms, selection: Collection[int], temperature: float
):
    # TODO importing a private function here... reimplement?
    reset = _maxwellboltzmanndistribution(
        atoms[selection].get_masses(), temperature * units.kB
    )
    _apply_momentum_to_selection(atoms, selection, reset)


def _scale_momentum_of_selection(
    atoms: Atoms, selection: Collection[int], temperature: float
):
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
