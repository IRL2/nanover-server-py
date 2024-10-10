import time
from unittest.mock import patch

import numpy
import pytest

from nanover.app import NanoverImdClient
from nanover.omni.omni import Simulation, OmniRunner
from nanover.testing import assert_equal_soon
from test_openmm import example_openmm
from test_ase_omm import example_ase_omm
from test_playback import example_playback
from openmm_simulation_utils import single_atom_simulation
from common import app_server

SIMULATION_FIXTURES = (
    "example_openmm",
    "example_playback",
    "example_ase_omm",
)

SIMULATION_FIXTURES_WITHOUT_PLAYBACK = [
    "example_openmm",
    "example_ase_omm",
]


TIMING_TOLERANCE = 0.05


@pytest.fixture
def multi_sim_runner(request):
    with OmniRunner.with_basic_server(port=0) as runner:
        for sim_fixture in SIMULATION_FIXTURES:
            simulation = request.getfixturevalue(sim_fixture)
            simulation.name = sim_fixture
            runner.add_simulation(simulation)
        yield runner


@pytest.fixture
def multi_sim_client_runner(multi_sim_runner):
    with NanoverImdClient.connect_to_single_server(
        port=multi_sim_runner.app_server.port
    ) as client:
        yield client, multi_sim_runner


@pytest.fixture
def multi_sim_client_runner_without_playback(request):
    with OmniRunner.with_basic_server(port=0) as runner:
        for sim_fixture in SIMULATION_FIXTURES_WITHOUT_PLAYBACK:
            simulation = request.getfixturevalue(sim_fixture)
            simulation.name = sim_fixture
            runner.add_simulation(simulation)
        with NanoverImdClient.connect_to_single_server(
            port=runner.app_server.port
        ) as client:
            yield client, runner


@pytest.mark.parametrize("sim_fixture", SIMULATION_FIXTURES)
def test_reset(sim_fixture, app_server, request):
    sim: Simulation = request.getfixturevalue(sim_fixture)
    sim.reset(app_server)


@pytest.mark.parametrize("sim_fixture", SIMULATION_FIXTURES)
def test_load(sim_fixture, request):
    sim: Simulation = request.getfixturevalue(sim_fixture)
    sim.load()


# TODO: not true for playback because stepping might step through state changes...
@pytest.mark.parametrize("sim_fixture", SIMULATION_FIXTURES_WITHOUT_PLAYBACK)
def test_step_gives_exactly_one_frame(sim_fixture, request, app_server):
    """
    Test that stepping sends exactly one frame.
    """
    sim: Simulation = request.getfixturevalue(sim_fixture)
    sim.load()
    sim.reset(app_server)
    sim.advance_by_one_step()

    with patch.object(
        app_server.frame_publisher, "send_frame", autospec=True
    ) as send_frame:
        for i in range(1, 20):
            sim.advance_by_one_step()
            assert send_frame.call_count == i


def test_list_simulations(multi_sim_runner):
    """
    Test runner lists exactly the expected simulations.
    """
    names = multi_sim_runner.list()["simulations"]
    assert sorted(names) == sorted(SIMULATION_FIXTURES)


def test_next_simulation_increments_counter(multi_sim_client_runner_without_playback):
    """
    Test each next command increments the simulation counter to the correct value.
    """
    client, runner = multi_sim_client_runner_without_playback
    client.subscribe_to_frames()

    for i in range(5):
        client.run_next()
        client.wait_until_first_frame()
        assert_equal_soon(lambda: client.current_frame.simulation_counter, lambda: i)


@pytest.mark.parametrize("fps", (5, 10, 30))
def test_play_step_interval(multi_sim_client_runner_without_playback, fps):
    """
    Test that the play step interval is respected and the sent frame frequency matches the requested interval.
    We only guarantee that the interval is close on average.
    """
    # We need at least a few frames to see intervals between
    test_frames = 30

    play_step_interval = 1 / fps
    client, runner = multi_sim_client_runner_without_playback

    client.subscribe_to_all_frames()
    runner.next()
    runner.runner.play_step_interval = play_step_interval

    while len(client.frames) < test_frames:
        time.sleep(0.1)
    runner.pause()

    # first frame (topology) isn't subject to intervals
    timestamps = [frame.server_timestamp for frame in client.frames[1:]]
    deltas = numpy.diff(timestamps)

    # The interval is not very accurate. We only check that the observed
    # interval is close on average.
    assert numpy.average(deltas) == pytest.approx(
        play_step_interval, abs=TIMING_TOLERANCE
    )
