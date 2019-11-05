import time

import pytest
from google.protobuf.struct_pb2 import Struct
from grpc import RpcError
from mock import Mock

from narupa.command.command_client import CommandClient
from narupa.command.command_server import CommandServer


@pytest.fixture
def client():
    with CommandClient() as c:
        yield c


@pytest.fixture
def server():
    with CommandServer() as s:
        yield s


def test_get_commands(client, server):
    mock = Mock()
    server.service.register_command("test", mock.callback, Struct())

    commands = client.get_commands_blocking()
    assert len(commands) == 1
    assert commands[0].name == "test"


def test_get_multiple_commands(client, server):
    mock = Mock()
    expected_names = set()
    for i in range(10):
        name = str(i)
        expected_names.add(name)
        server.service.register_command(name, mock.callback, Struct())

    commands = client.get_commands_blocking()
    assert len(commands) == 10
    assert set([command.name for command in commands]) == expected_names
    assert all(command.arguments == Struct() for command in commands)


def test_get_command_with_argument(client, server):
    mock = Mock()
    arguments = Struct()
    arguments['x'] = 1
    arguments['y'] = 2
    server.service.register_command("test", mock.callback, arguments)
    commands = client.get_commands_blocking()
    assert commands[0].arguments == arguments


def test_run_command(client, server):
    mock = Mock()
    server.service.register_command("test", mock.callback, Struct())
    client.run_command("test", Struct())
    time.sleep(0.1)
    assert mock.callback.call_count == 1


def test_run_no_args(client, server):
    mock = Mock()
    server.service.register_command("test", mock.callback)
    client.run_command("test")
    time.sleep(0.1)
    assert mock.callback.call_count == 1


def test_run_command_with_args(client, server):
    def sqr(struct: Struct) -> Struct:
        x = struct['x']
        result = Struct()
        result['y'] = x * x
        return result

    example_params = Struct()
    example_params['x'] = 2
    server.service.register_command("sqr", sqr, example_params)
    reply = client.run_command("sqr", example_params)
    assert reply.result['y'] == 4


def test_unknown_command(client, server):
    with pytest.raises(RpcError):
        client.run_command("unknown", Struct())
