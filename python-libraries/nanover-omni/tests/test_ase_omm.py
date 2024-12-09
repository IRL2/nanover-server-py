import warnings

import pytest
from ase import units

from nanover.ase.wall_constraint import VelocityWallConstraint
from nanover.omni.ase_omm import ASEOpenMMSimulation, CONSTRAINTS_UNSUPPORTED_MESSAGE

from common import app_server, ARGON_XML_PATH
from nanover.openmm.serializer import deserialize_simulation

from openmm_simulation_utils import (
    basic_system,
    basic_simulation,
    basic_simulation_with_imd_force,
    BASIC_SIMULATION_POSITIONS,
    empty_imd_force,
    assert_basic_simulation_topology,
    single_atom_system,
    single_atom_simulation,
    single_atom_simulation_with_imd_force,
    ARGON_SIMULATION_POSITION,
    assert_single_atom_simulation_topology,
)


@pytest.fixture
def example_ase_omm(app_server, single_atom_simulation):
    sim = ASEOpenMMSimulation.from_simulation(single_atom_simulation)
    sim.load()
    sim.reset(app_server)
    yield sim


def test_step_interval(example_ase_omm):
    """
    Test that advancing by one step increments the dynamics steps by frame_interval.
    """
    for i in range(5):
        assert (
            example_ase_omm.dynamics.get_number_of_steps()
            == i * example_ase_omm.frame_interval
        )
        example_ase_omm.advance_by_one_step()


# TODO: test it actually outputs
def test_verbose(example_ase_omm, app_server):
    """
    Test verbose option steps without exceptions.
    """
    example_ase_omm.verbose = True
    example_ase_omm.reset(app_server)
    for i in range(5):
        example_ase_omm.advance_by_one_step()


@pytest.mark.parametrize("walls", (False, True))
def test_walls(example_ase_omm, walls):
    example_ase_omm.use_walls = walls
    example_ase_omm.load()
    assert (
        any(
            isinstance(constraint, VelocityWallConstraint)
            for constraint in example_ase_omm.atoms.constraints
        )
        == walls
    )


def test_no_constraint_no_warning(example_ase_omm, app_server):
    """
    Test that a system without constraints does not cause a constraint warning
    to be logged.
    """
    with warnings.catch_warnings():
        warnings.simplefilter("error", UserWarning)
        example_ase_omm.load()
        example_ase_omm.reset(app_server)


def test_constraint_warning(example_ase_omm, app_server, recwarn):
    """
    Test that a system with constraints causes a constraint warning to be
    logged.
    """
    with pytest.warns(UserWarning, match=CONSTRAINTS_UNSUPPORTED_MESSAGE):
        example_ase_omm.simulation.system.addConstraint(0, 1, 1)
        example_ase_omm.load()
        example_ase_omm.reset(app_server)
