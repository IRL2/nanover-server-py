from unittest.mock import Mock

import pytest
from nanover.core.commands import CommandService


@pytest.fixture
def service():
    return CommandService(add_list_command=False)


def test_register_command(service):
    mock = Mock()
    service.register_command("test", mock.callback)
    assert len(service.commands) == 1
    assert "test" in service.commands


def test_register_command_with_args(service):
    mock = Mock()
    service.register_command("test", mock.callback, default_arguments={"a": 2})
    assert len(service.commands) == 1
    assert service.commands["test"].arguments == {"a": 2}


def test_unregister_command(service):
    mock = Mock()
    service.register_command("test", mock.callback)
    assert len(service.commands) == 1
    service.unregister_command("test")
    assert len(service.commands) == 0
    assert "test" not in service.commands


def test_unregister_nonexisting_command(service):
    service.unregister_command("test")
