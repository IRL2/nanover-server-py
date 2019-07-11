import pytest
from ase import units
from narupa.ase.openmm import OpenMMIMDRunner
from narupa.ase.openmm.runner import ImdParams
from narupa.trajectory.frame_server import DEFAULT_PORT as TRAJ_DEFAULT_PORT
from narupa.imd.imd_server import DEFAULT_PORT as IMD_DEFAULT_PORT

from .simulation_utils import basic_simulation, serialized_simulation_path


@pytest.fixture()
def runner(basic_simulation):
    runner = OpenMMIMDRunner(basic_simulation)
    yield runner
    runner.close()


def test_from_xml(serialized_simulation_path):
    with OpenMMIMDRunner.from_xml(serialized_simulation_path) as runner:
        assert runner.simulation is not None


def test_defaults(runner):
    default_params = ImdParams()
    assert runner.verbose == default_params.verbose
    assert runner.frame_interval == default_params.frame_interval
    assert runner.time_step == default_params.time_step
    assert runner.imd_port == default_params.imd_port
    assert runner.trajectory_port == default_params.trajectory_port
    assert runner.address == default_params.address


def test_dynamics_initialised(runner):
    assert runner.dynamics is not None


def test_run(runner):
    runner.run(10)
    assert runner.dynamics.get_number_of_steps() == 10


def test_frames_sent(runner):
    runner.run(12)
    assert runner.imd.frame_server.frame_count == 2


def test_verbose(basic_simulation):
    params = ImdParams(verbose=True)
    with OpenMMIMDRunner(basic_simulation, params) as runner:
        runner.run(10)


def test_frame_interval(basic_simulation):
    params = ImdParams(frame_interval=1)
    with OpenMMIMDRunner(basic_simulation, params) as runner:
        runner.run(2)
        assert runner.imd.frame_server.frame_count == 2


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
        runner = OpenMMIMDRunner(basic_simulation, params)
    # If the runner could be built, we need to close it. Though, at this point,
    # the runner may not exist.
    try:
        runner.close()
    except NameError:
        pass
