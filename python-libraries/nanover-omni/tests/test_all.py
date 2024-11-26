import time
from contextlib import contextmanager
from unittest.mock import patch

import numpy
import pytest

from nanover.app import NanoverImdClient
from nanover.omni.omni import Simulation, CLEAR_PREFIXES
from nanover.testing import assert_equal_soon, assert_in_soon, assert_not_in_soon
from nanover.utilities.change_buffers import DictionaryChange
from test_openmm import example_openmm
from test_ase_omm import example_ase_omm
from test_ase import example_ase, example_dynamics
from test_playback import example_playback
from openmm_simulation_utils import single_atom_simulation
from common import make_runner, make_connected_client_from_runner, make_app_server

SIMULATION_FIXTURES = (
    "example_openmm",
    "example_playback",
    "example_ase_omm",
    "example_ase",
)

SIMULATION_FIXTURES_WITHOUT_PLAYBACK = [
    "example_openmm",
    "example_ase_omm",
    "example_ase",
]


SIMULATION_FIXTURES_IMD = [
    "example_openmm",
    "example_ase_omm",
    "example_ase",
]


TIMING_TOLERANCE = 0.05


def make_sim_from_fixture(request, fixture_name):
    simulation: Simulation = request.getfixturevalue(fixture_name)
    simulation.name = fixture_name
    return simulation


def make_imd_sims(request):
    return [make_sim_from_fixture(request, fixture_name) for fixture_name in SIMULATION_FIXTURES_IMD]


def make_all_sims(request):
    return [make_sim_from_fixture(request, fixture_name) for fixture_name in SIMULATION_FIXTURES]


@pytest.fixture
def runner_with_all_sims(request):
    with make_runner(make_all_sims(request)) as runner:
        yield runner


@pytest.fixture
def runner_with_imd_sims(request):
    with make_runner(make_imd_sims(request)) as runner:
        yield runner


@pytest.mark.parametrize("sim_fixture", SIMULATION_FIXTURES)
def test_reset(sim_fixture, request):
    sim = make_sim_from_fixture(request, sim_fixture)

    with make_app_server() as app_server:
        sim.reset(app_server)


@pytest.mark.parametrize("sim_fixture", SIMULATION_FIXTURES)
def test_load(sim_fixture, request):
    sim = make_sim_from_fixture(request, sim_fixture)
    sim.load()


# TODO: not true for playback because stepping might step through state changes...
@pytest.mark.parametrize("sim_fixture", SIMULATION_FIXTURES_WITHOUT_PLAYBACK)
def test_step_gives_exactly_one_frame(sim_fixture, request):
    """
    Test that stepping sends exactly one frame.
    """
    with make_app_server() as app_server:
        sim = make_sim_from_fixture(request, sim_fixture)
        sim.load()
        sim.reset(app_server)
        sim.advance_by_one_step()

        with patch.object(
            app_server.frame_publisher, "send_frame", autospec=True
        ) as send_frame:
            for i in range(1, 20):
                sim.advance_by_one_step()
                assert send_frame.call_count == i


def test_list_simulations(runner_with_all_sims):
    """
    Test runner lists exactly the expected simulations.
    """
    names = runner_with_all_sims.list()["simulations"]
    assert sorted(names) == sorted(SIMULATION_FIXTURES)


def test_next_simulation_increments_counter(runner_with_imd_sims):
    """
    Test each next command increments the simulation counter to the correct value.
    """
    with make_connected_client_from_runner(runner_with_imd_sims) as client:
        client.subscribe_to_frames()

        for i in range(5):
            client.run_next()
            client.wait_until_first_frame()
            assert_equal_soon(lambda: client.current_frame.simulation_counter, lambda: i)


@pytest.mark.parametrize("fps", (5, 10, 30))
def test_play_step_interval(runner_with_imd_sims, fps):
    """
    Test that the play step interval is respected and the sent frame frequency matches the requested interval.
    We only guarantee that the interval is close on average.
    """
    # We need at least a few frames to see intervals between
    test_frames = 30

    play_step_interval = 1 / fps

    with make_connected_client_from_runner(runner_with_imd_sims) as client:
        client.subscribe_to_all_frames()
        runner_with_imd_sims.next()
        runner_with_imd_sims.runner.play_step_interval = play_step_interval

        while len(client.frames) < test_frames:
            time.sleep(0.1)
        runner_with_imd_sims.pause()

        # first frame (topology) isn't subject to intervals
        timestamps = [frame.server_timestamp for frame in client.frames[1:]]
        deltas = numpy.diff(timestamps)

    # The interval is not very accurate. We only check that the observed
    # interval is close on average.
    assert numpy.average(deltas) == pytest.approx(
        play_step_interval, abs=TIMING_TOLERANCE
    )


def test_simulation_switch_clears_state(runner_with_all_sims):
    """
    Test that state keys with certain prefixes are no longer present in the state after switching simulation.
    """
    key = "pytest"

    updates = {prefix + key: {} for prefix in CLEAR_PREFIXES}
    locks = {key: 10 for key in updates}

    with make_connected_client_from_runner(runner_with_all_sims) as client:
        client.attempt_update_multiplayer_state(DictionaryChange(updates=updates))
        client.attempt_update_multiplayer_locks(locks)

    with make_connected_client_from_runner(runner_with_all_sims) as client:
        client.subscribe_multiplayer()

        for key in updates:
            assert_in_soon(lambda: key, lambda: client._multiplayer_client.copy_state())

        client.run_next()

        for key in updates:
            assert_not_in_soon(
                lambda: key, lambda: client._multiplayer_client.copy_state()
            )
