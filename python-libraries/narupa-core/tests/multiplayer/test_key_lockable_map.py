from narupa.multiplayer.key_lockable_map import KeyLockableMap, ResourceLockedException
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
    with pytest.raises(ResourceLockedException):
        key_map.delete("2", key)
