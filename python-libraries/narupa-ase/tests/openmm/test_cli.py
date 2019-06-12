import pytest

from .simulation_utils import basic_simulation, serialized_simulation_path
from narupa.ase.openmm.cli import initialise


def test_initialise(serialized_simulation_path):
    args = [str(serialized_simulation_path)]
    runner = initialise(args)
    assert runner.simulation is not None


def test_timestep(serialized_simulation_path):
    args = [str(serialized_simulation_path), '-s', '0.5']
    runner = initialise(args)
    assert runner.time_step == 0.5


def test_interval(serialized_simulation_path):
    args = [str(serialized_simulation_path), '-f', '2']
    runner = initialise(args)
    assert runner.frame_interval == 2


def test_address(serialized_simulation_path):
    args = [str(serialized_simulation_path), '-a', 'localhost']
    runner = initialise(args)
    assert runner.address == 'localhost'


def test_traj_port(serialized_simulation_path):
    args = [str(serialized_simulation_path), '-t', '54324']
    runner = initialise(args)
    assert runner.trajectory_port == 54324


def test_imd_port(serialized_simulation_path):
    args = [str(serialized_simulation_path), '-i', '54324']
    runner = initialise(args)
    assert runner.imd_port == 54324


def test_same_port(serialized_simulation_path):
    args = [str(serialized_simulation_path), '-t', '54324', '-i', '54324']
    with pytest.raises(ValueError):
        _ = initialise(args)
