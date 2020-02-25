# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.
import os
from logging import Handler, WARNING

import pytest
from ase import units
from narupa.app.app_server import DEFAULT_NARUPA_PORT
from ase.io import read
from narupa.ase.openmm import OpenMMIMDRunner
from narupa.ase.openmm.runner import ImdParams, CONSTRAINTS_UNSUPPORTED_MESSAGE, LoggingParams
from narupa.core import DEFAULT_SERVE_ADDRESS
from narupa.essd import DiscoveryClient
from narupa.trajectory.frame_server import DEFAULT_PORT as TRAJ_DEFAULT_PORT
from narupa.imd.imd_server import DEFAULT_PORT as IMD_DEFAULT_PORT
from narupa.multiplayer.multiplayer_server import DEFAULT_PORT as MULTIPLAYER_DEFAULT_PORT
from narupa.ase.openmm.calculator import OpenMMCalculator
from narupa.ase.wall_calculator import VelocityWallCalculator

from .simulation_utils import basic_simulation, serialized_simulation_path

TRAJECTORY_OUTPUT_FILENAME = 'test.xyz'


class ListLogHandler(Handler):
    """
    Records records from a logger into a list that can be examined for testing.
    """

    def __init__(self):
        super().__init__()
        self.records = []

    def emit(self, record):
        self.records.append(record)

    def count_records(self, message: str, levelno: int):
        return sum(1 for record in self.records if record.message == message and record.levelno == levelno)


@pytest.fixture()
def params():
    """
    Test ImdParams set to use any available port, to avoid port clashes between tests.
    """
    params = ImdParams(
        port=0
    )
    return params


@pytest.fixture()
def logging_params(tmp_path):
    params = LoggingParams(
        trajectory_file=os.path.join(tmp_path, TRAJECTORY_OUTPUT_FILENAME),
        write_interval=1,
    )
    return params


@pytest.fixture()
def runner(basic_simulation, params):
    with OpenMMIMDRunner(basic_simulation, imd_params=params) as runner:
        yield runner


@pytest.fixture()
def default_runner(basic_simulation):
    with OpenMMIMDRunner(basic_simulation) as runner:
        yield runner


def test_from_xml(serialized_simulation_path, params):
    with OpenMMIMDRunner.from_xml(serialized_simulation_path,
                                  params=params) as runner:
        assert runner.simulation is not None


@pytest.mark.serial
def test_defaults(default_runner):
    runner = default_runner
    default_params = ImdParams()
    assert runner.verbose == default_params.verbose
    assert runner.frame_interval == default_params.frame_interval
    assert runner.time_step == default_params.time_step
    assert runner.port == DEFAULT_NARUPA_PORT
    assert runner.address == DEFAULT_SERVE_ADDRESS


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
    assert runner.app_server.frame_publisher.last_frame_index > 0


def test_verbose(basic_simulation, params):
    params.verbose = True
    with OpenMMIMDRunner(basic_simulation, params) as runner:
        runner.run(10)


@pytest.mark.parametrize('interval', (1, 2, 3))
def test_frame_interval(basic_simulation, interval, params):
    """
    Test that the frame server receives frames at the correct interval of
    dynamics steps.
    """
    params.frame_interval = interval
    with OpenMMIMDRunner(basic_simulation, params) as runner:
        runner.run(1)
        prev = runner.app_server.frame_publisher.last_frame_index
        runner.run(interval * 3)
        assert runner.app_server.frame_publisher.last_frame_index == prev + 3


def test_time_step(basic_simulation, params):
    params.time_step = 0.5
    with OpenMMIMDRunner(basic_simulation, params) as runner:
        assert runner.dynamics.dt == pytest.approx(0.5 * units.fs)


@pytest.mark.parametrize('walls, expected_calculator_class', (
        (False, OpenMMCalculator),
        (True, VelocityWallCalculator),
))
def test_walls(basic_simulation, walls, expected_calculator_class, params):
    params.walls = walls

    with OpenMMIMDRunner(basic_simulation, params) as runner:
        assert isinstance(runner._md_calculator, expected_calculator_class)


def test_no_constraint_no_warning(basic_simulation, params):
    """
    Test that a system without constraints does not cause a constraint warning
    to be logged.
    """
    handler = ListLogHandler()

    with OpenMMIMDRunner(basic_simulation, params) as runner:
        runner._logger.addHandler(handler)
        runner._validate_simulation()
        assert handler.count_records(CONSTRAINTS_UNSUPPORTED_MESSAGE, WARNING) == 0


def test_constraint_warning(basic_simulation, params):
    """
    Test that a system with constraints causes a constraint warning to be
    logged.
    """
    handler = ListLogHandler()
    basic_simulation.system.addConstraint(0, 1, 1)

    with OpenMMIMDRunner(basic_simulation, params) as runner:
        runner._logger.addHandler(handler)
        runner._validate_simulation()
        assert handler.count_records(CONSTRAINTS_UNSUPPORTED_MESSAGE, WARNING) == 1


def test_no_discovery(basic_simulation, params):
    params.discovery = False
    with OpenMMIMDRunner(basic_simulation, params) as runner:
        assert not runner.running_discovery


@pytest.mark.serial
def test_discovery(basic_simulation, params):
    with OpenMMIMDRunner(basic_simulation, params) as runner:
        assert runner.running_discovery
        assert runner.app_server.discovery is not None
        assert len(runner.app_server.discovery.services) == 1


@pytest.mark.serial
def test_discovery_with_client(basic_simulation, params):
    params.name = 'ASE Test Runner'
    with OpenMMIMDRunner(basic_simulation, params) as runner:
        assert runner.running_discovery
        discovery = runner.app_server.discovery
        assert len(discovery.services) == 1
        with DiscoveryClient() as client:
            # There may be servers running already, we only want to look at the
            # one we created in that test. We select it by name.
            servers = set(client.search_for_services(search_time=0.8, interval=0.01))
            relevant_servers = [server for server in servers
                                if server.name == params.name]
            assert len(relevant_servers) == 1
            server = relevant_servers[0]
            assert server in discovery.services
            assert server.name == runner.name
            assert len(server.services) == 3
            for port in server.services.values():
                assert port == runner.app_server.port


def test_logging(basic_simulation, params, logging_params):
    with OpenMMIMDRunner(basic_simulation, params, logging_params) as runner:
        runner.run(1)

    trajectory_file = runner.logging_info.trajectory_path
    assert os.path.exists(trajectory_file)
    assert runner.logging_info.base_trajectory_path == logging_params.trajectory_file

    atom_images = read(trajectory_file, index=':')
    assert len(atom_images) == 2


def test_no_logging(basic_simulation, params, logging_params):
    with OpenMMIMDRunner(basic_simulation, params) as runner:
        assert runner.logging_info is None


def test_logging_rate(basic_simulation, params, logging_params):
    logging_params.write_interval = 10
    expected_frames = 3

    with OpenMMIMDRunner(basic_simulation, params, logging_params) as runner:
        runner.run(expected_frames * logging_params.write_interval)

    trajectory_file = runner.logging_info.trajectory_path
    assert os.path.exists(trajectory_file)

    atom_images = read(trajectory_file, index=':')
    # always get one more frame as it dumps out the initial frame
    assert len(atom_images) == expected_frames + 1
