import warnings

import pytest

from nanover.ase.wall_constraint import VelocityWallConstraint
from nanover.omni.ase_omm import ASEOpenMMSimulation, CONSTRAINTS_UNSUPPORTED_MESSAGE

from common import make_app_server, make_loaded_sim

from openmm_simulation_utils import build_single_atom_simulation


@pytest.fixture
def ase_omm_single_atom():
    yield make_ase_omm_single_atom()


def make_ase_omm_single_atom():
    return ASEOpenMMSimulation.from_simulation(build_single_atom_simulation())


def test_step_interval(ase_omm_single_atom):
    """
    Test that advancing by one step increments the dynamics steps by frame_interval.
    """
    with make_loaded_sim(ase_omm_single_atom):
        for i in range(5):
            assert (
                ase_omm_single_atom.dynamics.get_number_of_steps()
                == i * ase_omm_single_atom.frame_interval
            )
            ase_omm_single_atom.advance_by_one_step()


# TODO: test it actually outputs
def test_verbose(ase_omm_single_atom):
    """
    Test verbose option steps without exceptions.
    """
    ase_omm_single_atom.verbose = True
    with make_loaded_sim(ase_omm_single_atom):
        for i in range(5):
            ase_omm_single_atom.advance_by_one_step()


@pytest.mark.parametrize("walls", (False, True))
def test_walls(ase_omm_single_atom, walls):
    ase_omm_single_atom.use_walls = walls
    ase_omm_single_atom.load()
    assert (
        any(
            isinstance(constraint, VelocityWallConstraint)
            for constraint in ase_omm_single_atom.atoms.constraints
        )
        == walls
    )


def test_no_constraint_no_warning(ase_omm_single_atom):
    """
    Test that a system without constraints does not cause a constraint warning
    to be logged.
    """
    with warnings.catch_warnings():
        warnings.simplefilter("error", UserWarning)
        with make_loaded_sim(ase_omm_single_atom):
            pass


def test_constraint_warning(ase_omm_single_atom, recwarn):
    """
    Test that a system with constraints causes a constraint warning to be
    logged.
    """
    with pytest.warns(UserWarning, match=CONSTRAINTS_UNSUPPORTED_MESSAGE):
        ase_omm_single_atom.simulation.system.addConstraint(0, 1, 1)
        with make_loaded_sim(ase_omm_single_atom):
            pass
