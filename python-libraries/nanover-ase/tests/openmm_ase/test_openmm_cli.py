import os

import pytest

from .simulation_utils import basic_simulation, serialized_simulation_path
from nanover.ase.openmm.cli import initialise
from nanover.ase.openmm.calculator import OpenMMCalculator
from nanover.ase.wall_constraint import VelocityWallConstraint


@pytest.fixture
def any_port():
    return ["-p", "0"]


def test_initialise(serialized_simulation_path, any_port):
    args = [str(serialized_simulation_path)] + any_port
    with initialise(args) as runner:
        assert runner.simulation is not None


def test_timestep(serialized_simulation_path, any_port):
    args = [str(serialized_simulation_path), "-s", "0.5"] + any_port
    with initialise(args) as runner:
        assert runner.simulation.time_step == 0.5


def test_interval(serialized_simulation_path, any_port):
    args = [str(serialized_simulation_path), "-f", "2"] + any_port
    with initialise(args) as runner:
        assert runner.simulation.frame_interval == 2


@pytest.mark.serial
def test_port(serialized_simulation_path):
    PORT = 29070  # The port reserved for Jedi Knight: Jedi Academy (2003), so should be safe.
    args = [str(serialized_simulation_path)] + ["-p", str(PORT)]
    with initialise(args) as runner:
        assert runner.app_server.port == PORT


def test_name(serialized_simulation_path, any_port):
    args = [str(serialized_simulation_path), "--name", "Test Server"] + any_port
    with initialise(args) as runner:
        assert runner.app_server.name == "Test Server"


@pytest.mark.parametrize(
    "wall_argument, has_walls",
    (
        ("-w", True),
        ("--walls", True),
        (None, False),
    ),
)
def test_walls(serialized_simulation_path, wall_argument, has_walls, any_port):
    args = [str(serialized_simulation_path)] + any_port
    if wall_argument is not None:
        args.append(wall_argument)
    with initialise(args) as runner:
        runner.simulation.load()
        assert (
            any(
                isinstance(constraint, VelocityWallConstraint)
                for constraint in runner.simulation.atoms.constraints
            )
            == has_walls
        )
