import time
from typing import Tuple

import pytest
from grpc import RpcError
from mock import Mock

from nanover.core.nanover_client import NanoverClient
from nanover.core.nanover_server import NanoverServer

TEST_COMMAND_KEY = "test"
MULTIPLY_KEY = "multiply"


@pytest.fixture
def client_server() -> Tuple[NanoverClient, NanoverServer]:
    with NanoverServer(address="localhost", port=0) as server:
        with NanoverClient.insecure_channel(
            address="localhost", port=server.port
        ) as client:
            yield client, server


@pytest.fixture
def default_args():
    return {"a": 2, "b": [1, 3, 4], "c": True}


@pytest.fixture
def mock_callback(default_args):
    return Mock(return_value=default_args)


def test_is_channel_owner(client_server):
    client, _ = client_server
    assert client.is_channel_owner


def test_not_channel_owner(client_server):
    client, _ = client_server
    with NanoverClient(channel=client.channel) as second_client:
        assert not second_client.is_channel_owner


def test_available_commands(client_server, default_args):
    """
    tests that the cached set of available commands matches those returned when updating.
    """
    client, server = client_server
    mock = Mock()
    server.register_command(TEST_COMMAND_KEY, mock.callback, default_args)
    commands = client.update_available_commands()

    assert client.available_commands == commands


def test_available_commands_init(client_server, default_args):
    """
    tests that the available commands is empty at initialisation of the client, if server has no
    commands.
    """
    client, server = client_server
    assert client.available_commands == {}


def test_client_inits_if_no_server():
    """
    tests that the client successfully initialises, even if there is no server.
    """
    with NanoverClient.insecure_channel(address="localhost", port=68393):
        pass


def test_get_commands(client_server, default_args):
    client, server = client_server
    mock = Mock()
    server.register_command(TEST_COMMAND_KEY, mock.callback, default_args)

    commands = client.update_available_commands()
    assert TEST_COMMAND_KEY in commands
    assert commands[TEST_COMMAND_KEY] == default_args


def test_commands_on_server(client_server):
    client, server = client_server
    mock = Mock()

    server.register_command(TEST_COMMAND_KEY, mock.callback)
    commands = server.commands
    assert TEST_COMMAND_KEY in commands
    command_registration = commands[TEST_COMMAND_KEY]
    assert command_registration.name == TEST_COMMAND_KEY
    assert command_registration.handler == mock.callback


def test_unregister_command(client_server):
    client, server = client_server
    mock = Mock()

    server.register_command(TEST_COMMAND_KEY, mock.callback)

    commands = client.update_available_commands()
    assert TEST_COMMAND_KEY in commands

    server.unregister_command(TEST_COMMAND_KEY)
    commands = client.update_available_commands()
    assert TEST_COMMAND_KEY not in commands


def test_get_multiple_commands(client_server):
    client, server = client_server
    mock = Mock()
    expected_names = set()
    for i in range(10):
        name = str(i)
        expected_names.add(name)
        server.register_command(name, mock.callback)

    commands = client.update_available_commands()
    assert set(commands.keys()).issuperset(expected_names)
    assert all(arguments == {} for arguments in commands.values())


def test_get_command_with_argument(client_server):
    client, server = client_server
    mock = Mock()
    arguments = {"x": 1, "y": 2}
    server.register_command(TEST_COMMAND_KEY, mock.callback, arguments)
    commands = client.update_available_commands()
    assert commands[TEST_COMMAND_KEY] == arguments


def test_run_command(client_server, mock_callback):
    client, server = client_server
    server.register_command(TEST_COMMAND_KEY, mock_callback)
    results = client.run_command(TEST_COMMAND_KEY, **{})
    time.sleep(0.1)
    mock_callback.assert_called_once()
    assert results == mock_callback()


def test_run_no_args(client_server, mock_callback):
    client, server = client_server
    server.register_command(TEST_COMMAND_KEY, mock_callback)
    client.run_command(TEST_COMMAND_KEY)
    time.sleep(0.1)
    mock_callback.assert_called_once()


def multiply(x=2, z=1):
    result = {"y": x * z}
    return result


def multiply_positional(x, z):
    result = {"y": x * z}
    return result


@pytest.mark.parametrize(
    "args, expected_result",
    [
        ({"x": 8}, 8),
        ({"z": 2}, 8),
        ({"x": 4, "z": 4}, 16),
        ({}, 4),
    ],
)
def test_run_command_with_args(client_server, args, expected_result):
    """
    tests that for various combinations of overriding default arguments, the expected result is returned.
    """
    client, server = client_server
    example_params = {"x": 4, "z": 1}
    server.register_command(MULTIPLY_KEY, multiply, example_params)
    reply = client.run_command(MULTIPLY_KEY, **args)
    assert {"y": expected_result} == reply


@pytest.mark.parametrize(
    "args, expected_result",
    [
        ({"x": 8}, 8),
        ({"z": 2}, 4),
        ({"x": 4, "z": 4}, 16),
        ({}, 2),
    ],
)
def test_run_command_with_args(client_server, args, expected_result):
    """
    tests that a command that takes arguments, but has no default args set upon registration, works correctly
    """
    client, server = client_server
    server.register_command(MULTIPLY_KEY, multiply)
    reply = client.run_command(MULTIPLY_KEY, **args)
    assert {"y": expected_result} == reply


@pytest.mark.parametrize(
    "args",
    [
        {"invalid_arg": 8},
        {"x": {"key": "value"}},
    ],
)
def test_run_command_with_invalid_args(client_server, args):
    """
    tests that running a command with invalid arguments raises exceptions.
    """
    client, server = client_server
    server.register_command(MULTIPLY_KEY, multiply)
    with pytest.raises(RpcError):
        _ = client.run_command(MULTIPLY_KEY, **args)


def test_run_command_with_positional_args(client_server):
    """
    tests that running a command with positional arguments raises an exception.
    """
    client, server = client_server
    server.register_command(MULTIPLY_KEY, multiply_positional)
    with pytest.raises(RpcError):
        _ = client.run_command(MULTIPLY_KEY)


def test_run_command_with_no_result(client_server):
    def method():
        pass

    client, server = client_server
    server.register_command(TEST_COMMAND_KEY, method)
    reply = client.run_command(TEST_COMMAND_KEY)
    assert reply == {}


def test_unknown_command(client_server):
    client, server = client_server
    with pytest.raises(RpcError):
        client.run_command("unknown")


def test_command_invalid_results(client_server):
    """
    tests that a command that has unserialisable results
    throws an rpc exception.
    """

    def invalid_method():
        result = {"y": object()}
        return result

    client, server = client_server
    server.register_command(TEST_COMMAND_KEY, invalid_method)
    with pytest.raises(RpcError):
        client.run_command(TEST_COMMAND_KEY)
