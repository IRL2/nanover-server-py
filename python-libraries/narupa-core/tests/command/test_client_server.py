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


def test_get_commands(client_server):
    client, server = client_server
    mock = Mock()
    server.register_command("test", mock.callback, Struct())

    commands = client.get_commands_blocking()
    assert len(commands) == 1
    assert commands[0].name == "test"


def test_get_multiple_commands(client_server):
    client, server = client_server
    mock = Mock()
    expected_names = set()
    for i in range(10):
        name = str(i)
        expected_names.add(name)
        server.register_command(name, mock.callback, Struct())

    commands = client.get_commands_blocking()
    assert len(commands) == 10
    assert set([command.name for command in commands]) == expected_names
    assert all(command.arguments == Struct() for command in commands)


def test_get_command_with_argument(client_server):
    client, server = client_server
    mock = Mock()
    arguments = Struct()
    arguments['x'] = 1
    arguments['y'] = 2
    server.register_command("test", mock.callback, arguments)
    commands = client.get_commands_blocking()
    assert commands[0].arguments == arguments


def test_run_command(client_server):
    client, server = client_server
    mock = Mock()
    server.register_command("test", mock.callback, Struct())
    client.run_command("test", Struct())
    time.sleep(0.1)
    assert mock.callback.call_count == 1


def test_run_no_args(client_server):
    client, server = client_server
    mock = Mock()
    server.register_command("test", mock.callback)
    client.run_command("test")
    time.sleep(0.1)
    assert mock.callback.call_count == 1


def test_run_command_with_args(client_server):
    def sqr(struct: Struct) -> Struct:
        x = struct['x']
        result = Struct()
        result['y'] = x * x
        return result

    client, server = client_server
    example_params = Struct()
    example_params['x'] = 2
    server.register_command("sqr", sqr, example_params)
    reply = client.run_command("sqr", example_params)
    assert reply.result['y'] == 4


def test_unknown_command(client_server):
    client, server = client_server
    with pytest.raises(RpcError):
        client.run_command("unknown", Struct())
