import pytest
from ase import units
from narupa.ase.openmm import OpenMMIMDRunner
from narupa.ase.openmm.runner import ImdParams
from narupa.trajectory.frame_server import DEFAULT_PORT as TRAJ_DEFAULT_PORT
from narupa.imd.imd_server import DEFAULT_PORT as IMD_DEFAULT_PORT
from narupa.ase.openmm.calculator import OpenMMCalculator
from narupa.ase.wall_calculator import VelocityWallCalculator

from .simulation_utils import basic_simulation, serialized_simulation_path


@pytest.fixture()
def runner(basic_simulation):
    with OpenMMIMDRunner(basic_simulation) as runner:
        yield runner


def test_from_xml(serialized_simulation_path):
    with OpenMMIMDRunner.from_xml(serialized_simulation_path) as runner:
        assert runner.simulation is not None


def test_defaults(runner):
    default_params = ImdParams()
    assert runner.verbose == default_params.verbose
    assert runner.frame_interval == default_params.frame_interval
    assert runner.time_step == default_params.time_step
    assert runner.imd_port == IMD_DEFAULT_PORT
    assert runner.trajectory_port == TRAJ_DEFAULT_PORT
    assert runner.address == default_params.address


def test_dynamics_initialised(runner):
    assert runner.dynamics is not None


def test_run(runner):
    runner.run(10)
    assert runner.dynamics.get_number_of_steps() == 10


def test_frames_sent(runner):
    """
    Test that the frame server has received frames after running dynamics.
    """
    runner.run(12)
    assert runner.imd.frame_server.frame_count > 0


def test_verbose(basic_simulation):
    params = ImdParams(verbose=True)
    with OpenMMIMDRunner(basic_simulation, params) as runner:
        runner.run(10)


@pytest.mark.parametrize('interval', (1, 2, 3))
def test_frame_interval(basic_simulation, interval):
    """
    Test that the frame server receives frames at the correct interval of
    dynamics steps.
    """
    params = ImdParams(frame_interval=interval)
    with OpenMMIMDRunner(basic_simulation, params) as runner:
        runner.run(1)
        prev = runner.imd.frame_server.frame_count
        runner.run(interval * 3)
        assert runner.imd.frame_server.frame_count == prev + 3


def test_time_step(basic_simulation):
    params = ImdParams(time_step=0.5)
    with OpenMMIMDRunner(basic_simulation, params) as runner:
        assert runner.dynamics.get_timestep() == pytest.approx(0.5 * units.fs)


@pytest.mark.parametrize('trajectory_port, imd_port', (
    (5555, 5555),
    (None, TRAJ_DEFAULT_PORT),
    (IMD_DEFAULT_PORT, None),
))
def test_same_port_failure(basic_simulation, trajectory_port, imd_port):
    params = ImdParams()
    params.address = 'localhost'
    params.trajectory_port = trajectory_port
    params.imd_port = imd_port

    with pytest.raises(ValueError):
        with OpenMMIMDRunner(basic_simulation, params):
            pass


@pytest.mark.parametrize('trajectory_port, imd_port', (
        (0, 0),
        (0, None),
        (None, 0),
        (0, 80),
        (80, 0),
))
def test_same_port_accepts_zero(trajectory_port, imd_port):
    assert not OpenMMIMDRunner._services_use_same_port(trajectory_port, imd_port)


@pytest.mark.parametrize('walls, expected_calculator_class', (
    (False, OpenMMCalculator),
    (True, VelocityWallCalculator),
))
def test_walls(basic_simulation, walls, expected_calculator_class):
    params = ImdParams()
    params.trajectory_port = 0
    params.imd_port = 0
    params.walls = walls

    with OpenMMIMDRunner(basic_simulation, params) as runner:
        assert isinstance(runner._md_calculator, expected_calculator_class)
