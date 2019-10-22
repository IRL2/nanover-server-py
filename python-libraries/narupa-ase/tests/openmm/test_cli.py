# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.
import pytest

from .simulation_utils import basic_simulation, serialized_simulation_path
from narupa.ase.openmm.cli import initialise
from narupa.ase.openmm.calculator import OpenMMCalculator
from narupa.ase.wall_calculator import VelocityWallCalculator


def test_initialise(serialized_simulation_path):
    args = [str(serialized_simulation_path)]
    with initialise(args) as runner:
        assert runner.simulation is not None


def test_timestep(serialized_simulation_path):
    args = [str(serialized_simulation_path), '-s', '0.5']
    with initialise(args) as runner:
        assert runner.time_step == 0.5


def test_interval(serialized_simulation_path):
    args = [str(serialized_simulation_path), '-f', '2']
    with initialise(args) as runner:
        assert runner.frame_interval == 2


def test_address(serialized_simulation_path):
    args = [str(serialized_simulation_path), '-a', 'localhost']
    with initialise(args) as runner:
        assert runner.address == 'localhost'


def test_traj_port(serialized_simulation_path):
    args = [str(serialized_simulation_path), '-t', '54324']
    with initialise(args) as runner:
        assert runner.trajectory_port == 54324


def test_imd_port(serialized_simulation_path):
    args = [str(serialized_simulation_path), '-i', '54324']
    with initialise(args) as runner:
        assert runner.imd_port == 54324


def test_discovery_service(serialized_simulation_path):
    args = [str(serialized_simulation_path)]
    with initialise(args) as runner:
        assert runner.running_discovery is True


def test_discovery_service_not_running(serialized_simulation_path):
    args = [str(serialized_simulation_path), '--no_discovery']
    with initialise(args) as runner:
        assert runner.running_discovery is False


def test_discovery_service_port(serialized_simulation_path):
    args = [str(serialized_simulation_path), '--discovery_port', '88888']
    with initialise(args) as runner:
        assert runner.discovery_port == 88888


def test_name(serialized_simulation_path):
    args = [str(serialized_simulation_path), '--name', 'Test Server']
    with initialise(args) as runner:
        assert runner.name == 'Test Server'


def test_discovery_service_port_not_running(serialized_simulation_path):
    args = [str(serialized_simulation_path), '--no_discovery']
    with initialise(args) as runner:
        with pytest.raises(AttributeError):
            _ = runner.discovery_port


def test_same_port(serialized_simulation_path):
    args = [str(serialized_simulation_path), '-t', '54324', '-i', '54324']
    with pytest.raises(ValueError):
        _ = initialise(args)


@pytest.mark.parametrize('wall_argument, expected_calculator_class', (
        ('-w', VelocityWallCalculator),
        ('--walls', VelocityWallCalculator),
        (None, OpenMMCalculator),
))
def test_walls(serialized_simulation_path, wall_argument, expected_calculator_class):
    args = [str(serialized_simulation_path)]
    if wall_argument is not None:
        args.append(wall_argument)
    with initialise(args) as runner:
        assert isinstance(runner._md_calculator, expected_calculator_class)
