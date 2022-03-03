# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.
import contextlib
from typing import ContextManager, Tuple

import pytest
from ase import units
from ase.lattice.bravais import Lattice
from ase.lattice.cubic import FaceCenteredCubic
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.calculators.lj import LennardJones
from ase.md import Langevin
from ase.md.md import MolecularDynamics
from ase.md.nvtberendsen import NVTBerendsen
from hypothesis import strategies, given
from narupa.ase.imd_calculator import (ImdCalculator, _get_cancelled_interactions, _get_atoms_to_reset,
                                       _scale_momentum_of_selection)
from narupa.imd import ImdServer
from narupa.imd.particle_interaction import ParticleInteraction
from narupa.trajectory.frame_data import MissingDataError
import numpy as np
from test_imd_calculator import imd_calculator_co
from util import imd_server

# Sets of dummy interactions to test cancellation selections.
NUM_INTERACTIONS = 10

INTERACTIONS_NO_RESET = {
    ("0", str(i)):
        ParticleInteraction(particles=(10 * i, 10 * i + 1))
    for i in range(NUM_INTERACTIONS)}
INTERACTIONS_RESET = {
    (str(i), "1"):
        ParticleInteraction(reset_velocities=True, particles=(i, i + 1))
    for i in range(1, NUM_INTERACTIONS + 1)}

# the set of atoms that should be reset based on atoms selected in INTERACTIONS_RESET
ATOMS_TO_RESET = set(i for i in range(1, NUM_INTERACTIONS + 2))
ALL_INTERACTIONS = dict(INTERACTIONS_NO_RESET)
ALL_INTERACTIONS.update(INTERACTIONS_RESET)

MIN_TEMP = 1e-6
MAX_TEMP = 1e10
TEST_TEMPERATURE = 300


def fcc_atoms():
    """
    A FCC crystal for use as a set of atoms in tests.
    :return: ASE atoms object.
    """
    size = 2
    # Set up a crystal
    atoms = FaceCenteredCubic(directions=[[1, 0, 0], [0, 1, 0], [0, 0, 1]],
                              symbol="Cu",
                              size=(size, size, size),
                              pbc=True)
    return atoms


@contextlib.contextmanager
def imd_calculator_berendsen_dynamics_context() -> ContextManager[Tuple[ImdCalculator, Lattice, NVTBerendsen]]:
    server = ImdServer(address=None, port=0)
    atoms = fcc_atoms()
    calculator = LennardJones()
    dynamics = NVTBerendsen(atoms, 1.0, TEST_TEMPERATURE, 1.0)
    imd_calculator = ImdCalculator(server.imd_state, calculator, atoms, dynamics=dynamics)
    yield imd_calculator, atoms, dynamics
    server.close()


@pytest.fixture
def imd_calculator_berendsen_dynamics():
    """
    Initialises an IMD calculator with berendsen NVT integrator and an FCC crystal.
    """
    with imd_calculator_berendsen_dynamics_context() as imd_calculator:
        yield imd_calculator


@pytest.fixture
def imd_calculator_langevin_dynamics():
    """
    Initialises an IMD calculator with langevin NVT integrator and an FCC crystal.
    """
    server = ImdServer(address=None, port=0)
    atoms = fcc_atoms()
    calculator = LennardJones()
    dynamics = Langevin(atoms, 1.0, friction=1.0, temperature_K=TEST_TEMPERATURE)
    imd_calculator = ImdCalculator(server.imd_state, calculator, atoms, dynamics=dynamics)
    yield imd_calculator, atoms, dynamics
    server.close()


@pytest.mark.parametrize("new, old, exp_cancelled",
                         [(INTERACTIONS_NO_RESET, INTERACTIONS_RESET, INTERACTIONS_RESET),
                          (INTERACTIONS_NO_RESET, INTERACTIONS_NO_RESET, {}),
                          ({}, INTERACTIONS_NO_RESET, INTERACTIONS_NO_RESET),
                          (ALL_INTERACTIONS, INTERACTIONS_NO_RESET, {}),
                          (INTERACTIONS_RESET, ALL_INTERACTIONS, INTERACTIONS_NO_RESET)
                          ]
                         )
def test_cancelled_interactions(new, old, exp_cancelled):
    """
    Tests that the expected set of interaction keys are produced when cancelling interactions.
    """
    cancelled = _get_cancelled_interactions(new, old)
    assert cancelled.keys() == exp_cancelled.keys()


@pytest.mark.parametrize("interactions, exp_atoms",
                         [(INTERACTIONS_NO_RESET, set()),
                          (INTERACTIONS_RESET, ATOMS_TO_RESET),
                          (ALL_INTERACTIONS, ATOMS_TO_RESET),
                          ({}, set()),
                          ])
def test_get_atoms_to_reset(interactions, exp_atoms):
    """
    Tests that the expected set of atoms to reset are produced when passing cancelled interactions.
    A set of cancelled interactions may include both those marked to reset atoms and those not,
    but only the atoms to reset are returned.
    """
    atoms = _get_atoms_to_reset(interactions)
    assert atoms == exp_atoms


def test_temperature_not_set(imd_calculator_co):
    """
    Tests handling of temperature not set in an IMD calculator.
    """
    calculator, atoms, _ = imd_calculator_co
    with pytest.raises(MissingDataError):
        _ = calculator.temperature


def test_custom_temperature():
    """
    Tests handling temperature not set in constructor of IMD calculator.
    """
    server = ImdServer(address=None, port=0)
    atoms = fcc_atoms()
    calculator = LennardJones()
    imd_calculator = ImdCalculator(server.imd_state, calculator, atoms, reset_scale=0.1)
    imd_calculator.temperature = 100
    assert pytest.approx(imd_calculator.reset_temperature) == 0.1 * 100


def test_temperature_not_set_md(imd_calculator_co):
    """
    Tests handling of temperature not set in an IMD calculator with a dynamics
    object that does not implement temperature.
    """
    calculator, atoms, _ = imd_calculator_co
    # molecular dynamics object does not implement a temperature by default.
    calculator._dynamics = MolecularDynamics(atoms, 1.0, None)
    with pytest.raises(MissingDataError):
        _ = calculator.temperature


def test_temperature_berendsen(imd_calculator_berendsen_dynamics):
    """
    Tests that berendsen NVT dynamics produces a temperature that
    can be used by IMD.
    """
    calculator, atoms, dynamics = imd_calculator_berendsen_dynamics
    assert calculator.temperature == dynamics.temperature


def test_temperature_langevin(imd_calculator_langevin_dynamics):
    """
    Tests that langevin NVT dynamics produces a temperature that
    can be used by IMD.
    """
    calculator, atoms, dynamics = imd_calculator_langevin_dynamics
    assert calculator.temperature == dynamics.temp


def test_temperature_custom(imd_calculator_co):
    """
    Tests that setting a custom temperature enables use of a temperature
    in IMD.
    """
    calculator, atoms, _ = imd_calculator_co
    calculator.temperature = TEST_TEMPERATURE
    assert calculator.temperature == TEST_TEMPERATURE


@strategies.composite
def random_atom_selection(draw):
    """
    From the FCC atoms, produces a random subset of indices.
    :return: Random subset of indices of FCC atoms, and the FCC atoms.
    """
    atoms = fcc_atoms()
    selection = strategies.sets(strategies.integers(0, len(atoms) - 1), min_size=1)
    return draw(selection), atoms


@given(random_atom_selection(), strategies.floats(MIN_TEMP, MAX_TEMP), strategies.floats(MIN_TEMP, MAX_TEMP))
def test_temperature_scaling_selection(random_atom_selection, dynamics_temp, reset_temperature):
    """
    Tests that scaling momentum in an atoms object on a selection of atoms produces the expected
    temperature on the selection, and leaves the other atoms untouched.
    """

    atoms_to_reset, atoms = random_atom_selection
    atoms_to_reset = np.array(list(atoms_to_reset))

    MaxwellBoltzmannDistribution(atoms, temperature_K=dynamics_temp, force_temp=True)
    assert pytest.approx(atoms.get_temperature()) == dynamics_temp

    # get the velocities of the unselected atoms before resetting.
    not_selected = inverse_selection(atoms, atoms_to_reset)
    velocities = atoms[not_selected].get_velocities()

    # reset atoms.
    _scale_momentum_of_selection(atoms, atoms_to_reset, reset_temperature)

    assert pytest.approx(atoms[atoms_to_reset].get_temperature()) == reset_temperature
    assert np.allclose(velocities, atoms[not_selected].get_velocities())


def test_no_reset(imd_calculator_berendsen_dynamics):
    calculator, atoms, dyn = imd_calculator_berendsen_dynamics
    MaxwellBoltzmannDistribution(atoms, temperature_K=300)
    velocities = atoms.get_velocities()
    calculator._reset_velocities(atoms, {}, INTERACTIONS_NO_RESET)
    reset_velocities = atoms.get_velocities()
    assert np.allclose(velocities, reset_velocities)


def generate_interactions(selection, num_interactions=1):
    particles_per_interaction = int(len(selection) / num_interactions)
    interactions = {}
    for i in range(num_interactions):
        particles = [selection[idx] for idx in
                     range(i * particles_per_interaction, (i + 1) * particles_per_interaction)]
        interaction = ParticleInteraction(particles=particles, reset_velocities=True)
        interactions[("0", str(i))] = interaction
    return interactions


def inverse_selection(collection, selection):
    mask = np.ones(len(collection), bool)
    mask[selection] = 0
    return mask


@given(atom_selection=random_atom_selection())
def test_reset_velocities(atom_selection):
    selection, _ = atom_selection
    selection = np.array(list(selection))
    with imd_calculator_berendsen_dynamics_context() as imd_calculator:
        calculator, atoms, dyn = imd_calculator
        MaxwellBoltzmannDistribution(atoms, temperature_K=300)

        not_selected = inverse_selection(atoms, selection)
        velocities = atoms[not_selected].get_velocities()

        interactions = generate_interactions(selection)
        calculator._reset_velocities(atoms, {}, interactions)
        assert pytest.approx(atoms[selection].get_temperature()) == calculator.reset_temperature
        assert np.allclose(velocities, atoms[not_selected].get_velocities())


def test_reset_calculator(imd_calculator_berendsen_dynamics):
    """
    Integration test checking that calling calculate twice, with
    an interaction then without, produces a reset in velocities.
    """
    imd_calculator = imd_calculator_berendsen_dynamics

    calculator, atoms, dyn = imd_calculator
    atoms.calc = calculator

    calculator.atoms = None
    calculator.calculate(atoms=atoms, system_changes=['numbers'])
    MaxwellBoltzmannDistribution(atoms, temperature_K=300)

    selection = [0, 1]
    selection = np.array(list(selection))
    interaction = ParticleInteraction(
        particles=selection,
        reset_velocities=True,
    )
    calculator._imd_state.insert_interaction('interaction.test', interaction)
    atoms.get_forces()
    calculator._imd_state.remove_interaction('interaction.test')
    atoms.get_forces()

    assert pytest.approx(atoms[selection].get_temperature()) == calculator.reset_temperature
