# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.
import os

import pytest

from .simulation_utils import basic_simulation, serialized_simulation_path
from narupa.ase.openmm.cli import initialise
from narupa.ase.openmm.calculator import OpenMMCalculator
from narupa.ase.wall_constraint import VelocityWallConstraint


@pytest.fixture
def any_port():
    return ['-p', '0']


def test_initialise(serialized_simulation_path, any_port):
    args = [str(serialized_simulation_path)] + any_port
    with initialise(args) as runner:
        assert runner.simulation is not None


def test_timestep(serialized_simulation_path, any_port):
    args = [str(serialized_simulation_path), '-s', '0.5'] + any_port
    with initialise(args) as runner:
        assert runner.time_step == 0.5


def test_interval(serialized_simulation_path, any_port):
    args = [str(serialized_simulation_path), '-f', '2'] + any_port
    with initialise(args) as runner:
        assert runner.frame_interval == 2


def test_address(serialized_simulation_path, any_port):
    # cannot run discovery here, as the CI servers cannot broadcast on localhost
    args = [str(serialized_simulation_path), '-a', 'localhost', '--no-discovery'] + any_port
    with initialise(args) as runner:
        assert runner.address == 'localhost'


@pytest.mark.serial
def test_port(serialized_simulation_path):
    PORT = 29070  # The port reserved for Jedi Knight: Jedi Academy (2003), so should be safe.
    args = [str(serialized_simulation_path)] + ['-p', str(PORT)]
    with initialise(args) as runner:
        assert runner.port == PORT


def test_discovery_service(serialized_simulation_path, any_port):
    args = [str(serialized_simulation_path)] + any_port
    with initialise(args) as runner:
        assert runner.running_discovery is True


def test_discovery_service_not_running(serialized_simulation_path, any_port):
    args = [str(serialized_simulation_path), '--no-discovery'] + any_port
    with initialise(args) as runner:
        assert runner.running_discovery is False


def test_discovery_service_port(serialized_simulation_path, any_port):
    args = [str(serialized_simulation_path), '--discovery-port', '88888'] + any_port
    with initialise(args) as runner:
        assert runner.discovery_port == 88888


def test_discovery_service_port_not_running(serialized_simulation_path, any_port):
    args = [str(serialized_simulation_path), '--no-discovery'] + any_port
    with initialise(args) as runner:
        with pytest.raises(AttributeError):
            _ = runner.discovery_port


def test_name(serialized_simulation_path, any_port):
    args = [str(serialized_simulation_path), '--name', 'Test Server'] + any_port
    with initialise(args) as runner:
        assert runner.name == 'Test Server'


@pytest.mark.parametrize('wall_argument, has_walls', (
        ('-w', True),
        ('--walls', True),
        (None, False),
))
def test_walls(serialized_simulation_path, wall_argument, has_walls, any_port):
    args = [str(serialized_simulation_path)] + any_port
    if wall_argument is not None:
        args.append(wall_argument)
    with initialise(args) as runner:
        assert any(isinstance(constraint, VelocityWallConstraint) for constraint in runner.atoms.constraints) == has_walls


@pytest.fixture()
def log_path(tmp_path):
    log_path = os.path.join(tmp_path, "test.xyz")
    return log_path


def test_trajectory_logging(serialized_simulation_path, log_path, any_port):
    args = [str(serialized_simulation_path), '-o', log_path] + any_port
    with initialise(args) as runner:
        assert runner.logging_info
        assert not os.path.exists(runner.logging_info.trajectory_path)
        runner.run(1)
        assert os.path.exists(runner.logging_info.trajectory_path)


def test_trajectory_no_logging(serialized_simulation_path, log_path, any_port):
    args = [str(serialized_simulation_path)] + any_port
    with initialise(args) as runner:
        assert runner.logging_info is None


def test_trajectory_logging_rate(serialized_simulation_path, log_path, any_port):
    args = [str(serialized_simulation_path), '-o', log_path, '--write-interval', '10'] + any_port
    with initialise(args) as runner:
        assert runner.logging_info.write_interval == 10
