# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.
import pytest

from .simulation_utils import basic_simulation, serialized_simulation_path
from narupa.ase.openmm.cli import initialise
from narupa.ase.openmm.calculator import OpenMMCalculator
from narupa.ase.wall_calculator import VelocityWallCalculator


@pytest.fixture
def test_ports():
    return ['-t', '0', '-i', '0', '-m', '0']


def set_port_arg(option_list, flag, port):
    flag_index = option_list.index(flag)
    value_index = flag_index + 1
    option_list[value_index] = port


def test_initialise(serialized_simulation_path, test_ports):
    args = [str(serialized_simulation_path)] + test_ports
    with initialise(args) as runner:
        assert runner.simulation is not None


def test_timestep(serialized_simulation_path, test_ports):
    args = [str(serialized_simulation_path), '-s', '0.5'] + test_ports
    with initialise(args) as runner:
        assert runner.time_step == 0.5


def test_interval(serialized_simulation_path, test_ports):
    args = [str(serialized_simulation_path), '-f', '2'] + test_ports
    with initialise(args) as runner:
        assert runner.frame_interval == 2


def test_address(serialized_simulation_path, test_ports):
    # cannot run discovery here, as the CI servers cannot broadcast on localhost
    args = [str(serialized_simulation_path), '-a', 'localhost', '--no-discovery'] + test_ports
    with initialise(args) as runner:
        assert runner.address == 'localhost'


def test_traj_port(serialized_simulation_path, test_ports):
    set_port_arg(test_ports, '-t', '62034')
    args = [str(serialized_simulation_path)] + test_ports
    with initialise(args) as runner:
        assert runner.trajectory_port == 62034


def test_imd_port(serialized_simulation_path, test_ports):
    set_port_arg(test_ports, '-i', '62035')
    args = [str(serialized_simulation_path)] + test_ports
    with initialise(args) as runner:
        assert runner.imd_port == 62035


def test_discovery_service(serialized_simulation_path, test_ports):
    args = [str(serialized_simulation_path)] + test_ports
    with initialise(args) as runner:
        assert runner.running_discovery is True


def test_discovery_service_not_running(serialized_simulation_path, test_ports):
    args = [str(serialized_simulation_path), '--no-discovery'] + test_ports
    with initialise(args) as runner:
        assert runner.running_discovery is False


def test_discovery_service_port(serialized_simulation_path, test_ports):
    args = [str(serialized_simulation_path), '--discovery-port', '88888'] + test_ports
    with initialise(args) as runner:
        assert runner.discovery_port == 88888


def test_discovery_service_port_not_running(serialized_simulation_path, test_ports):
    args = [str(serialized_simulation_path), '--no-discovery'] + test_ports
    with initialise(args) as runner:
        with pytest.raises(AttributeError):
            _ = runner.discovery_port


def test_name(serialized_simulation_path, test_ports):
    args = [str(serialized_simulation_path), '--name', 'Test Server'] + test_ports
    with initialise(args) as runner:
        assert runner.name == 'Test Server'


def test_multiplayer(serialized_simulation_path, test_ports):
    args = [str(serialized_simulation_path)] + test_ports
    with initialise(args) as runner:
        assert runner.running_multiplayer is True


def test_multiplayer_not_running(serialized_simulation_path, test_ports):
    args = [str(serialized_simulation_path), '--no-multiplayer'] + test_ports
    with initialise(args) as runner:
        assert runner.running_multiplayer is False


def test_multiplayer_port(serialized_simulation_path, test_ports):
    set_port_arg(test_ports, '-m', '64321')
    args = [str(serialized_simulation_path)] + test_ports
    with initialise(args) as runner:
        assert runner.multiplayer_port == 64321


def test_multiplayer_port_not_set(serialized_simulation_path, test_ports):
    args = [str(serialized_simulation_path), '--no-multiplayer'] + test_ports
    with initialise(args) as runner:
        with pytest.raises(AttributeError):
            _ = runner.multiplayer_port


def test_same_port(serialized_simulation_path, test_ports):
    set_port_arg(test_ports, '-t', '54324')
    set_port_arg(test_ports, '-i', '54324')
    args = [str(serialized_simulation_path)] + test_ports
    with pytest.raises(ValueError):
        _ = initialise(args)


@pytest.mark.parametrize('wall_argument, expected_calculator_class', (
        ('-w', VelocityWallCalculator),
        ('--walls', VelocityWallCalculator),
        (None, OpenMMCalculator),
))
def test_walls(serialized_simulation_path, wall_argument, expected_calculator_class, test_ports):
    args = [str(serialized_simulation_path)] + test_ports
    if wall_argument is not None:
        args.append(wall_argument)
    with initialise(args) as runner:
        assert isinstance(runner._md_calculator, expected_calculator_class)
