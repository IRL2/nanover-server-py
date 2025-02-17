import time
import pytest
import numpy as np
from ase.calculators.lj import LennardJones
from nanover.ase import converter
from nanover.ase.imd_calculator import (
    ImdCalculator,
    get_periodic_box_lengths,
)
from nanover.ase.null_calculator import NullCalculator
from nanover.imd import ImdClient
from nanover.imd.particle_interaction import ParticleInteraction
from util import co_atoms, imd_server, client_interaction, state_wrapper, c_atoms


@pytest.fixture
def atoms():
    return co_atoms()


@pytest.fixture
def interact_c():
    interaction = ParticleInteraction(
        position=[1, 0, 0],
        particles=[0],
        scale=100.0,
        interaction_type="spring",
    )
    return interaction


@pytest.fixture
def imd_calculator_co(imd_server):
    atoms = co_atoms()
    calculator = LennardJones()
    imd_calculator = ImdCalculator(imd_server.imd_state, calculator, atoms)
    yield imd_calculator, atoms, imd_server


def test_imd_calculator_no_interactions(imd_calculator_co):
    imd_calculator, atoms, _ = imd_calculator_co
    properties = ("energy", "forces")
    imd_calculator.calculator.calculate(atoms=atoms, properties=properties)
    expected_results = imd_calculator.calculator.results
    imd_calculator.calculate()
    results = imd_calculator.results
    for key in properties:
        assert np.allclose(results[key], expected_results[key])
    assert results["interactive_energy"] == 0
    assert np.all(results["interactive_forces"] == np.zeros((len(atoms), 3)))


def test_imd_calculator_shape_change_error(imd_calculator_co):
    imd_calculator, _, _ = imd_calculator_co

    with pytest.raises(AssertionError):
        atoms_different = c_atoms()
        imd_calculator.calculate(atoms=atoms_different)


def test_imd_calculator_one_dimension_pbc(state_wrapper):
    calculator = NullCalculator()
    atoms = co_atoms()
    atoms.set_pbc((True, False, False))
    with pytest.raises(NotImplementedError):
        ImdCalculator(state_wrapper, calculator, atoms)


def test_imd_calculator_no_pbc(imd_calculator_co):
    imd_calculator, atoms, _ = imd_calculator_co
    atoms.set_pbc((False, False, False))
    assert get_periodic_box_lengths(atoms) is None


def test_imd_calculator_not_orthorhombic(state_wrapper):
    calculator = NullCalculator()
    atoms = co_atoms()
    atoms.set_cell([1, 1, 1, 45, 45, 45])
    with pytest.raises(NotImplementedError):
        ImdCalculator(state_wrapper, calculator, atoms)


@pytest.mark.parametrize(
    "position, imd_energy, imd_forces",
    [
        ([1, 0, 0], 1, [2, 0, 0, 0, 0, 0]),
        ([3, 0, 0], 1, [-2, 0, 0, 0, 0, 0]),
        ([-1, 0, 0], 1, [-2, 0, 0, 0, 0, 0]),
        ([0, 1, 0], 1, [0, 2, 0, 0, 0, 0]),
        ([0, 3, 0], 1, [0, -2, 0, 0, 0, 0]),
        ([0, -1, 0], 1, [0, -2, 0, 0, 0, 0]),
        ([0, 0, 1], 1, [0, 0, 2, 0, 0, 0]),
        ([0, 0, 3], 1, [0, 0, -2, 0, 0, 0]),
        ([0, 0, -1], 1, [0, 0, -2, 0, 0, 0]),
        ([5, 0, 0], 1, [2, 0, 0, 0, 0, 0]),
    ],
)
def test_one_interaction(
    position, imd_energy, imd_forces, imd_calculator_co, interact_c
):
    """
    tests an interaction in several different positions, including periodic boundary positions.
    """
    imd_calculator, atoms, imd_server = imd_calculator_co
    # reshape the expected imd forces.
    imd_forces = np.array(imd_forces)
    imd_forces = imd_forces.reshape((2, 3))

    # perform the internal energy calculation.
    properties = ("energy", "forces", "interactive_energy", "interactive_forces")
    imd_calculator.calculator.calculate(atoms=atoms, properties=properties)
    internal_energy = imd_calculator.calculator.results["energy"]
    internal_forces = imd_calculator.calculator.results["forces"]

    # perform the calculation with interaction applied.
    interact_c.position = position
    with ImdClient.insecure_channel(port=imd_server.port) as imd_client:
        with client_interaction(imd_client, interact_c):
            time.sleep(0.1)
            assert len(imd_calculator.interactions) == 1
            # Update interactions
            imd_calculator.update_interactions()
            imd_calculator.calculate(properties=properties)
            results = imd_calculator.results

    # set up the expected energy and forces.
    expected_imd_energy_kjmol = interact_c.scale * imd_energy * atoms.get_masses()[0]
    expected_imd_energy = expected_imd_energy_kjmol * converter.KJMOL_TO_EV
    expected_imd_forces = (
        interact_c.scale
        * atoms.get_masses()[0]
        * (imd_forces * converter.KJMOL_TO_EV / converter.NM_TO_ANG)
    )
    expected_forces = internal_forces + expected_imd_forces
    expected_energy = internal_energy + expected_imd_energy

    forces = results["forces"]
    assert np.allclose(results["interactive_energy"], expected_imd_energy)
    assert np.allclose(results["interactive_forces"], expected_imd_forces)
    assert np.allclose(results["energy"], expected_energy)
    assert np.allclose(forces, expected_forces)
