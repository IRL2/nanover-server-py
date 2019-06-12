import pytest
from ase import units
from narupa.ase.openmm import OpenMMIMDRunner
from narupa.ase.openmm.runner import ImdParams

from .simulation_utils import basic_simulation, serialized_simulation_path


def test_from_xml(serialized_simulation_path):
    runner = OpenMMIMDRunner.from_xml(serialized_simulation_path)
    assert runner.simulation is not None
    runner.close()


def test_defaults(basic_simulation):
    runner = OpenMMIMDRunner(basic_simulation)
    default_params = ImdParams()
    assert runner.verbose == default_params.verbose
    assert runner.frame_interval == default_params.frame_interval
    assert runner.time_step == default_params.time_step
    assert runner.imd_port == default_params.imd_port
    assert runner.trajectory_port == default_params.trajectory_port
    assert runner.address == default_params.address
    runner.close()


def test_dynamics_initialised(basic_simulation):
    runner = OpenMMIMDRunner(basic_simulation)
    assert runner.dynamics is not None
    runner.close()


def test_run(basic_simulation):
    runner = OpenMMIMDRunner(basic_simulation)
    runner.run(10)
    assert runner.dynamics.get_number_of_steps() == 10
    runner.close()


def test_frames_sent(basic_simulation):
    runner = OpenMMIMDRunner(basic_simulation)
    runner.run(12)
    assert runner.imd.frame_server.frame_count == 2
    runner.close()


def test_verbose(basic_simulation):
    params = ImdParams(verbose=True)
    runner = OpenMMIMDRunner(basic_simulation, params)
    runner.run(10)
    runner.close()


def test_frame_interval(basic_simulation):
    params = ImdParams(frame_interval=1)
    runner = OpenMMIMDRunner(basic_simulation, params)
    runner.run(2)
    assert runner.imd.frame_server.frame_count == 2
    runner.close()


def test_time_step(basic_simulation):
    params = ImdParams(time_step=0.5)
    runner = OpenMMIMDRunner(basic_simulation, params)
    assert runner.dynamics.get_timestep() == pytest.approx(0.5 * units.fs)
    runner.close()
