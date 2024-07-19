from unittest.mock import patch

import pytest

from nanover.app import NanoverImdClient
from nanover.omni.omni import Simulation, OmniRunner
from nanover.testing import assert_equal_soon
from test_openmm import example_openmm
from test_ase_omm import example_ase_omm
from test_playback import example_playback
from common import app_server

SIMULATION_FIXTURES = (
    "example_openmm",
    "example_playback",
    "example_ase_omm",
)


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


@pytest.mark.parametrize("sim_fixture", SIMULATION_FIXTURES)
def test_reset(sim_fixture, app_server, request):
    sim: Simulation = request.getfixturevalue(sim_fixture)
    sim.reset(app_server)


@pytest.mark.parametrize("sim_fixture", SIMULATION_FIXTURES)
def test_load(sim_fixture, request):
    sim: Simulation = request.getfixturevalue(sim_fixture)
    sim.load()


# TODO: not true for playback because stepping might step through state changes...
@pytest.mark.parametrize("sim_fixture", set(SIMULATION_FIXTURES) - {"example_playback"})
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


def test_next_simulation_increments_counter(multi_sim_client_runner):
    """
    Test each next command increments the simulation counter to the correct value.
    """
    client, runner = multi_sim_client_runner
    client.subscribe_to_frames()

    for i in range(5):
        client.run_next()
        client.wait_until_first_frame()
        assert_equal_soon(lambda: client.current_frame.simulation_counter, lambda: i)


# TODO: actually support this lol
@pytest.xfail
@pytest.mark.serial
@pytest.mark.parametrize("fps", (5, 10, 30))
def test_throttling(client_runner, fps):
    """
    The runner uses the requested MD throttling.

    Here we make sure the runner throttles the dynamics according to the
    dynamics interval. However, we only guarantee that the target dynamics
    interval is close on average.
    """
    # We need at least a few frames to see intervals between
    test_frames = 30

    dynamics_interval = 1 / fps
    client, runner = client_runner
    runner.dynamics_interval = dynamics_interval

    client.subscribe_to_all_frames()
    runner.run()

    while len(client.frames) < test_frames:
        time.sleep(0.1)

    runner.imd.cancel_run(wait=True)

    # first frame (topology) isn't subject to intervals
    timestamps = [frame.server_timestamp for frame in client.frames[1:]]
    deltas = numpy.diff(timestamps)

    # The interval is not very accurate. We only check that the observed
    # interval is close on average.
    assert numpy.average(deltas) == pytest.approx(
        dynamics_interval, abs=TIMING_TOLERANCE
    )
