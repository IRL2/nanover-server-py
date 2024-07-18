import pytest
from ase import units

from nanover.ase.wall_constraint import VelocityWallConstraint
from nanover.omni.ase_omm import ASEOpenMMSimulation

from common import app_server, ARGON_XML_PATH


@pytest.fixture
def example_ase_omm(app_server):
    sim = ASEOpenMMSimulation.from_xml_path(ARGON_XML_PATH)
    sim.load()
    sim.reset(app_server)
    yield sim


def test_step_interval(example_ase_omm):
    """
    Test that advancing by one step increments the dynamics steps by frame_interval.
    """
    for i in range(5):
        assert example_ase_omm.dynamics.get_number_of_steps() == i * example_ase_omm.frame_interval
        example_ase_omm.advance_by_one_step()


@pytest.mark.parametrize("time_step", (0.5, 1.0, 1.5))
def test_time_step(example_ase_omm, time_step, app_server):
    example_ase_omm.time_step = time_step
    example_ase_omm.frame_interval = 1
    example_ase_omm.reset(app_server)
    for i in range(5):
        assert example_ase_omm.dynamics.dt == pytest.approx(time_step * units.fs)
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
    example_ase_omm.simulation = None
    example_ase_omm.use_walls = walls
    example_ase_omm.load()
    assert (
        any(
            isinstance(constraint, VelocityWallConstraint)
            for constraint in example_ase_omm.atoms.constraints
        )
        == walls
    )


@pytest.mark.xfail
def test_no_constraint_no_warning(example_ase_omm):
    """
    Test that a system without constraints does not cause a constraint warning
    to be logged.
    """
    handler = ListLogHandler()

    with ASEOpenMMRunner(basic_simulation, imd_params) as runner:
        runner._logger.addHandler(handler)
        runner._validate_simulation()
        assert handler.count_records(CONSTRAINTS_UNSUPPORTED_MESSAGE, WARNING) == 0


@pytest.mark.xfail
def test_constraint_warning(example_ase_omm):
    """
    Test that a system with constraints causes a constraint warning to be
    logged.
    """
    handler = ListLogHandler()
    basic_simulation.system.addConstraint(0, 1, 1)

    with ASEOpenMMRunner(basic_simulation, imd_params) as runner:
        runner._logger.addHandler(handler)
        runner._validate_simulation()
        assert handler.count_records(CONSTRAINTS_UNSUPPORTED_MESSAGE, WARNING) == 1
