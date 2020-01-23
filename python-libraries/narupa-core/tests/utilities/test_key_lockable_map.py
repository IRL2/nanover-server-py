from narupa.utilities.key_lockable_map import (
    KeyLockableMap, ResourceLockedError,
)
import pytest


@pytest.fixture
def key_map():
    return KeyLockableMap()


def test_delete_key(key_map):
    key = "name"
    key_map.set("1", key, 2)
    assert key_map.get(key) == 2
    key_map.delete("1", key)
    assert key_map.get(key) is None


def test_delete_key_locked(key_map):
    key = "name"
    key_map.set("1", key, 2)
    key_map.lock_key("1", key)
    with pytest.raises(ResourceLockedError):
        key_map.delete("2", key)


def test_delete_missing_key(key_map):
    with pytest.raises(KeyError):
        key_map.delete("1", "non-existant-key")


def test_set_no_replace(key_map):
    key = "name"
    key_map.set_no_replace(key, 2)
    assert key_map.get(key) == 2


def test_set_no_replace_attempt_replace(key_map):
    key = "name"
    key_map.set_no_replace(key, 2)
    assert key_map.get(key) == 2
    with pytest.raises(KeyError):
        key_map.set_no_replace(key, 3)


def test_set_attempt_replace(key_map):
    key = "name"
    key_map.set("1", key, 2)
    assert key_map.get(key) == 2
    with pytest.raises(KeyError):
        key_map.set_no_replace(key, 3)


def test_set_no_replace_then_set(key_map):
    key = "name"
    key_map.set_no_replace(key, 2)
    assert key_map.get(key) == 2
    key_map.set("1", key, 3)
    assert key_map.get(key) == 3
