import pytest
from mock import Mock

from narupa.command import CommandService


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
