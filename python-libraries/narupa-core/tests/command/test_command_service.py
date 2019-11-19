import pytest
from google.protobuf.internal.well_known_types import Struct
from google.protobuf.struct_pb2 import ListValue
from mock import Mock
from narupa.core.command_info import dict_to_struct

from narupa.core.command_service import CommandService


@pytest.fixture
def service():
    return CommandService()


def test_register_command(service):
    mock = Mock()
    service.register_command("test", mock.callback)
    assert len(service.commands) == 1
    assert "test" in service.commands


def test_register_command_with_args(service):
    mock = Mock()
    service.register_command("test", mock.callback, {'a': 2})
    assert len(service.commands) == 1
    assert service.commands['test'].info.arguments == {'a': 2}


def test_register_existing_command(service):
    mock = Mock()
    service.register_command("test", mock.callback)
    with pytest.raises(ValueError):
        service.register_command("test", mock.callback)


def test_unregister_command(service):
    mock = Mock()
    service.register_command("test", mock.callback)
    assert len(service.commands) == 1
    service.unregister_command("test")
    assert len(service.commands) == 0
    assert "test" not in service.commands


def test_unregister_nonexisting_command(service):
    with pytest.raises(KeyError):
        service.unregister_command("test")


def assert_struct_dictionary_equal(struct, dictionary):
    assert len(dictionary) == len(struct)
    for key in dictionary:
        if isinstance(struct[key], ListValue):
            assert pytest.approx(struct[key]) == dictionary[key]
        elif isinstance(struct[key], Struct):
            assert_struct_dictionary_equal(struct[key], dictionary[key])
        else:
            assert struct[key] == dictionary[key]


@pytest.mark.parametrize('dictionary',
                         [
                             {'int': 1, 'bool': True, 'list': [1.0, 2.4, 3.8], 'str': 'hello'},
                             {},
                             {'recursive': {'int': 5}}
                         ])
def test_dict_to_struct(dictionary):
    struct = dict_to_struct(dictionary)
    assert_struct_dictionary_equal(struct, dictionary)

@pytest.mark.parametrize('dictionary',
                         [
                             {object(): 1},
                             {1: object()},
                         ])
def test_dict_to_struct_invalid(dictionary):
    with pytest.raises(ValueError):
        _ = dict_to_struct(dictionary)