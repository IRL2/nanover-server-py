from unittest.mock import patch, call

import pytest

from nanover.omni.playback import PlaybackSimulation, SCENE_POSE_IDENTITY

from common import RECORDING_PATH_TRAJ, RECORDING_PATH_STATE, make_loaded_sim


@pytest.fixture
def example_playback():
    with make_loaded_sim(make_example_playback()) as sim:
        yield sim


def make_example_playback():
    return PlaybackSimulation(
        "nanotube-example-recording",
        traj=RECORDING_PATH_TRAJ,
        state=RECORDING_PATH_STATE,
    )


def test_step_gives_exactly_one_emit(example_playback):
    with patch.object(example_playback, "emit", autospec=True) as emit:
        for i in range(1, 100):
            example_playback.advance_by_one_step()
            assert emit.call_count == i


def test_step_loops(example_playback):
    """
    Test that stepping past the last entry loops back to the beginning, emitting the first entry as expected.
    """
    first_time, first_frame, first_update = example_playback.entries[0]
    example_playback.next_entry_index = len(example_playback.entries) - 1

    with patch.object(example_playback, "emit", autospec=True) as emit:
        example_playback.advance_by_one_step()
        emit.reset_mock()
        example_playback.advance_by_one_step()

    assert example_playback.next_entry_index == 1
    emit.assert_called_with(frame=first_frame, update=first_update)


def test_advance_seconds_emits_intermediate(example_playback):
    """
    Test that advancing by a certain duration emits exactly all intermediate entries.
    """
    with patch.object(example_playback, "emit", autospec=True) as emit:
        for _ in range(5):
            entries = example_playback.entries[
                example_playback.next_entry_index : example_playback.next_entry_index
                + 5
            ]
            calls = [call(frame=frame, update=update) for _, frame, update in entries]

            emit.reset_mock()
            last_time, _, _ = entries[-1]
            example_playback.advance_by_seconds(last_time - example_playback.time)

            assert emit.call_count == len(calls)
            emit.assert_has_calls(calls)


def test_playback_finds_changed_keys(example_playback):
    """
    Test that the keys chosen for reset are exactly those changed in the recording.
    """
    assert example_playback.changed_keys == {
        "avatar.2c85cd8f-580e-4d41-98e5-b6f71be2575b",
        "playarea.2c85cd8f-580e-4d41-98e5-b6f71be2575b",
        "update.index.2c85cd8f-580e-4d41-98e5-b6f71be2575b",
        "interaction.8615bf5a-763a-48c5-812f-2ca5a06589c1",
    }


def test_playback_reset_clears_keys(example_playback):
    """
    Test that reset clears the relevant keys.
    """

    count = len(example_playback.entries)
    half = count // 2

    # playback to the middle of the simulation where there is an avatar and an interaction
    for _ in range(half):
        example_playback.advance_to_next_entry()

    example_playback.reset(example_playback.app_server)

    # check all keys that should be removed are not present
    with example_playback.app_server.server.lock_state() as state:
        assert all(key not in state for key in example_playback.changed_keys)


def test_playback_starts_with_scene_identity(example_playback):
    """
    Test that reset results in scene identity.
    """
    with example_playback.app_server.server.lock_state() as state:
        assert state["scene"] == SCENE_POSE_IDENTITY
