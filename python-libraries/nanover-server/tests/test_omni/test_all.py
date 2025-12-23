from unittest.mock import patch

import numpy
import pytest
from nanover.imd import ParticleInteraction

from nanover.app.omni import CLEAR_PREFIXES
from nanover.testing import (
    assert_equal_soon,
    assert_in_soon,
    assert_not_in_soon,
)
from nanover.testing.asserts import assert_true_soon
from nanover.trajectory.keys import SIMULATION_EXCEPTION, SERVER_TIMESTAMP
from nanover.utilities.change_buffers import DictionaryChange
from .test_openmm import make_example_openmm
from .test_ase import make_example_ase
from .test_playback import make_example_playback
from .test_playback_mda import make_example_playback_mda
from .common import (
    make_runner,
    make_connected_client_from_runner,
    make_app_server,
)

SIMULATION_FACTORIES_IMD = [
    make_example_openmm,
    make_example_ase,
]

SIMULATION_FACTORIES_PLAYBACK = [
    make_example_playback,
    make_example_playback_mda,
]

SIMULATION_FACTORIES_ALL = SIMULATION_FACTORIES_IMD + SIMULATION_FACTORIES_PLAYBACK


TIMING_TOLERANCE = 0.05


def make_imd_sims():
    return [sim_factory() for sim_factory in SIMULATION_FACTORIES_IMD]


def make_all_sims():
    return [sim_factory() for sim_factory in SIMULATION_FACTORIES_ALL]


@pytest.fixture
def runner_with_all_sims():
    with make_runner(*make_all_sims()) as runner:
        yield runner


@pytest.fixture
def runner_with_imd_sims():
    with make_runner(*make_imd_sims()) as runner:
        yield runner


@pytest.mark.parametrize("sim_factory", SIMULATION_FACTORIES_ALL)
def test_load(sim_factory):
    sim = sim_factory()
    sim.load()


@pytest.mark.parametrize("sim_factory", SIMULATION_FACTORIES_ALL)
def test_reset(sim_factory):
    sim = sim_factory()

    with make_app_server() as app_server:
        sim.load()
        sim.reset(app_server)


# TODO: not true for playback because stepping might step through state changes...
@pytest.mark.parametrize("sim_factory", SIMULATION_FACTORIES_IMD)
def test_step_gives_exactly_one_frame(sim_factory):
    """
    Test that stepping sends exactly one frame.
    """
    with make_app_server() as app_server:
        sim = sim_factory()
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
    assert sorted(names) == sorted(sim.name for sim in make_all_sims())


def test_next_simulation_increments_counter(runner_with_imd_sims):
    """
    Test each next command increments the simulation counter to the correct value.
    """
    with make_connected_client_from_runner(runner_with_imd_sims) as client:
        for i in range(5):
            client.run_next()
            assert_equal_soon(
                lambda: client.current_frame.simulation_counter, lambda: i
            )


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
        runner_with_imd_sims.next()
        runner_with_imd_sims.runner.play_step_interval = play_step_interval

        assert_true_soon(
            lambda: len(client.frames) >= test_frames,
            timeout=play_step_interval * test_frames * 2,
            interval=play_step_interval * 0.5,
        )
        runner_with_imd_sims.pause()

        # first frame (topology) isn't subject to intervals
        timestamps = [frame[SERVER_TIMESTAMP] for frame in client.frames[1:]]
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

    with make_connected_client_from_runner(runner_with_all_sims) as client:
        client.update_state(DictionaryChange(updates=updates))

        with make_connected_client_from_runner(runner_with_all_sims) as client:
            for key in updates:
                assert_in_soon(
                    lambda: key,
                    lambda: client._state_dictionary.copy_content(),
                )

            client.run_next()

            for key in updates:
                assert_not_in_soon(
                    lambda: key,
                    lambda: client._state_dictionary.copy_content(),
                )


@pytest.mark.parametrize("sim_factory", SIMULATION_FACTORIES_IMD)
def test_wait_until_minimum_usable_frame(sim_factory):
    """
    Test that the minimum usable frame contains topology and position information.
    """
    with make_runner(sim_factory()) as runner:
        with make_connected_client_from_runner(runner) as client:
            runner.load(0)
            first_frame = client.wait_until_minimum_usable_frame()

            assert (
                len(first_frame.particle_positions) > 0
                and len(first_frame.particle_elements) > 0
            )


@pytest.mark.parametrize("sim_factory", SIMULATION_FACTORIES_IMD)
def test_interaction_invalid_particle_index(sim_factory):
    """
    Test that attempting to calculate iMD forces for interactions with out of bounds particles adds an exception to the
    frame and removes all interactions.
    """
    with make_runner(sim_factory()) as runner:
        with make_connected_client_from_runner(runner) as client:
            runner.load(0)
            frame = client.wait_until_minimum_usable_frame()

            interaction_id = client.start_interaction(
                ParticleInteraction(particles=[frame.particle_count + 10])
            )

            # exception exists in frame
            assert_in_soon(lambda: SIMULATION_EXCEPTION, lambda: client.current_frame)
            # interaction no longer exists in state
            assert_not_in_soon(
                lambda: interaction_id,
                lambda: runner.app_server.imd.active_interactions,
            )
