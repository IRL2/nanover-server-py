import time

import pytest
from google.protobuf.struct_pb2 import Struct
from grpc import RpcError
from mock import Mock

from narupa.core.narupa_client import NarupaClient
from narupa.core.narupa_server import NarupaServer


@pytest.fixture
def client_server():
    with NarupaServer(address="localhost", port=0) as s:
        with NarupaClient(address="localhost", port=s.port) as c:
            yield c, s


@pytest.fixture
def default_args():
    return {'a': 2, 'b': [1, 3, 4], 'c': True}


@pytest.fixture
def mock_callback(default_args):
    return Mock(return_value=default_args)


def test_get_commands(client_server, default_args):
    client, server = client_server
    mock = Mock()
    server.register_command("test", mock.callback, default_args)

    commands = client.update_available_commands()
    assert len(commands) == 1
    assert commands[0].name == "test"
    assert commands[0].arguments == default_args


def test_get_multiple_commands(client_server):
    client, server = client_server
    mock = Mock()
    expected_names = set()
    for i in range(10):
        name = str(i)
        expected_names.add(name)
        server.register_command(name, mock.callback, Struct())

    commands = client.update_available_commands()
    assert len(commands) == 10
    assert set([command.name for command in commands]) == expected_names
    assert all(command.arguments == {} for command in commands)


def test_get_command_with_argument(client_server):
    client, server = client_server
    mock = Mock()
    arguments = {'x': 1, 'y': 2}
    server.register_command("test", mock.callback, arguments)
    commands = client.update_available_commands()
    assert commands[0].arguments == arguments


def test_run_command(client_server, mock_callback):
    client, server = client_server
    server.register_command("test", mock_callback)
    results = client.run_command("test", **{})
    time.sleep(0.1)
    mock_callback.assert_called_once()
    assert results == mock_callback()


def test_run_no_args(client_server, mock_callback):
    client, server = client_server
    server.register_command("test", mock_callback)
    client.run_command("test")
    time.sleep(0.1)
    mock_callback.assert_called_once()


def test_run_command_with_args(client_server):
    def sqr(x=2):
        result = {'y': x * x}
        return result

    client, server = client_server
    example_params = {'x': 2}
    server.register_command("sqr", sqr, example_params)
    reply = client.run_command("sqr", **example_params)
    assert reply == {'y': 4}


def test_run_command_with_no_result(client_server):
    def method():
        pass

    client, server = client_server
    server.register_command("test", method)
    reply = client.run_command("test")
    assert reply == {}


def test_unknown_command(client_server):
    client, server = client_server
    with pytest.raises(RpcError):
        client.run_command("unknown")
