import pytest
from hypothesis import given, strategies as st
from mock import Mock

from nanover.testing.servers import (
    make_connected_server_client_setup,
    connect_client_to_server,
)
from nanover.testing.utilities import simplify_numpy
from nanover.utilities.change_buffers import DictionaryChange

from nanover.testing import assert_equal_soon, assert_not_in_soon, assert_in_soon
from nanover.testing.strategies import (
    command_arguments,
    state_updates,
    frames,
    user_ids,
)
from nanover.websocket.client.app_client import NanoverImdClient


def echo(**arguments):
    return arguments


@pytest.fixture(scope="module")
def reusable_setup():
    with make_connected_server_client_setup() as setup:
        setup.app.register_command("test/identity", echo)
        yield setup


@pytest.fixture(scope="module")
def reusable_setup_two_clients():
    with make_connected_server_client_setup() as setup:
        setup.client.register_command("test/client_echo", echo)
        with connect_client_to_server(setup.server) as client:
            yield setup, client


@given(frame=frames())
def test_websocket_sends_frame(reusable_setup, frame):
    reusable_setup.server_publish_frame(frame)
    reusable_setup.assert_frames_synced_soon()


@given(frame=frames())
def test_websocket_sends_frame_two_clients(reusable_setup_two_clients, frame):
    reusable_setup, client2 = reusable_setup_two_clients

    reusable_setup.server_publish_frame(frame)
    reusable_setup.assert_frames_synced_soon()

    assert_equal_soon(
        lambda: simplify_numpy(reusable_setup.client.current_frame.frame_dict),
        lambda: simplify_numpy(client2.current_frame.frame_dict),
    )

    reusable_setup.server_publish_frame_reset()
    reusable_setup.assert_frames_synced_soon()

    assert_equal_soon(
        lambda: simplify_numpy(reusable_setup.client.current_frame.frame_dict),
        lambda: simplify_numpy(client2.current_frame.frame_dict),
    )


@given(updates=state_updates())
def test_websocket_sends_state(reusable_setup, updates):
    """
    Test that state updates made directly on the server are accurately reflected on the client.
    """
    change = DictionaryChange(updates=updates)
    reusable_setup.server_update_state(change)

    assert_equal_soon(
        lambda: pick(reusable_setup.server_current_state, updates),
        lambda: updates,
    )

    reusable_setup.assert_states_synced_soon()


@given(arguments=command_arguments())
def test_echo_command(reusable_setup, arguments):
    """
    Test that preprepared identity command successfully returns unaltered and arbitrary arguments the command
    is called with.
    """
    mock_callback = Mock()
    reusable_setup.client.run_command("test/identity", arguments, mock_callback)

    assert_equal_soon(
        lambda: mock_callback.call_args and mock_callback.call_args.args[0],
        lambda: arguments,
    )


@given(arguments=command_arguments())
def test_websocket_register_command(reusable_setup_two_clients, arguments):
    """
    Test one client can call a function registered by another.
    """
    setup, extra_client = reusable_setup_two_clients

    assert (
        extra_client.run_command_blocking("test/client_echo", **arguments) == arguments
    )


@given(frame=frames())
def test_client_frame_reset(reusable_setup, frame):
    reusable_setup.server_publish_frame(frame)
    reusable_setup.assert_frames_synced_soon()

    reusable_setup.server_publish_frame_reset()
    reusable_setup.assert_frames_synced_soon()


@given(
    user_id=user_ids(),
    fields=st.sets(st.text(min_size=1, max_size=16), min_size=1, max_size=4),
)
def test_disconnect_cleans_owned_keys(reusable_setup, user_id, fields):
    """
    Test that client disconnecting causes server to clear keys "owned" by that client.
    """
    updates = {f"{field}.{user_id}": True for field in fields}
    updates[f"avatar.{user_id}"] = True  # field used to establish user_id

    first = next(iter(updates))

    with NanoverImdClient.from_app_server(reusable_setup.app) as client:
        client.update_state(DictionaryChange(updates=updates))

        assert_in_soon(lambda: first, lambda: reusable_setup.server_current_state)
        assert all(key in reusable_setup.server_current_state for key in updates)

    assert_not_in_soon(lambda: first, lambda: reusable_setup.server_current_state)
    assert all(key not in reusable_setup.server_current_state for key in updates)


def pick(dictionary, keys):
    return {key: value for key, value in dictionary.items() if key in keys}
