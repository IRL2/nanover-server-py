# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.

from time import sleep

import pytest

from narupa.multiplayer.multiplayer_lock import MultiplayerObjectLock


@pytest.fixture
def test_object():
    return object()


def test_lock(test_object):
    multiplayer_lock = MultiplayerObjectLock(test_object)
    player_id = "1"
    result = multiplayer_lock.try_lock(player_id)
    assert result == True
    assert multiplayer_lock.is_locked() == True


def test_lock_object(test_object):
    multiplayer_lock = MultiplayerObjectLock(test_object)
    assert multiplayer_lock.locked_object is test_object


def test_lock_ownership(test_object):
    multiplayer_lock = MultiplayerObjectLock(test_object)
    player_id = "1"
    multiplayer_lock.try_lock(player_id)
    result = multiplayer_lock.is_lock_owner(player_id)
    assert result == True


def test_get_lock_owner(test_object):
    multiplayer_lock = MultiplayerObjectLock(test_object)
    player_id = "1"
    multiplayer_lock.try_lock(player_id)
    result = multiplayer_lock.get_lock_owner()
    assert result == "1"


def test_lock_ownership_unlocked(test_object):
    multiplayer_lock = MultiplayerObjectLock(test_object)
    result = multiplayer_lock.is_lock_owner("1")
    result = multiplayer_lock.get_lock_owner()
    assert result is None


def test_lock_ownership_not_owner(test_object):
    multiplayer_lock = MultiplayerObjectLock(test_object)
    player_id = "1"
    multiplayer_lock.try_lock(player_id)
    result = multiplayer_lock.is_lock_owner("2")
    assert result == False


def test_unlock_unlocked(test_object):
    multiplayer_lock = MultiplayerObjectLock(test_object)
    player_id = "1"
    result = multiplayer_lock.try_unlock(player_id)
    assert result == True


def test_unlock(test_object):
    multiplayer_lock = MultiplayerObjectLock(test_object)
    player_id = "1"
    multiplayer_lock.try_lock(player_id)
    assert multiplayer_lock.is_locked()
    result = multiplayer_lock.try_unlock(player_id)
    assert result == True
    assert multiplayer_lock.is_locked() == False


def test_second_lock_attempt(test_object):
    multiplayer_lock = MultiplayerObjectLock(test_object)
    player_id = "1"
    multiplayer_lock.try_lock(player_id)
    assert multiplayer_lock.is_locked()

    player_id = "2"
    result = multiplayer_lock.try_lock(player_id)
    assert result == False


def test_second_unlock_attempt(test_object):
    multiplayer_lock = MultiplayerObjectLock(test_object)
    player_id = "1"
    multiplayer_lock.try_lock(player_id)
    assert multiplayer_lock.is_locked()

    player_id = "2"
    result = multiplayer_lock.try_unlock(player_id)
    assert result == False


def test_elapsed_lock(test_object):
    multiplayer_lock = MultiplayerObjectLock(test_object, max_lock_time=1)
    player_id = "1"
    multiplayer_lock.try_lock(player_id)
    assert multiplayer_lock.is_locked()

    sleep(1)

    player_id = "2"
    result = multiplayer_lock.try_lock(player_id)
    assert result == True


def test_elapsed_unlock(test_object):
    multiplayer_lock = MultiplayerObjectLock(test_object, max_lock_time=0.01)
    player_id = "1"
    multiplayer_lock.try_lock(player_id)
    assert multiplayer_lock.is_locked()

    sleep(0.02)

    player_id = "2"
    result = multiplayer_lock.try_unlock(player_id)
    assert result == True


def test_update_object(test_object):
    multiplayer_lock = MultiplayerObjectLock(test_object)
    player_id = "1"
    new_object = object()
    success = multiplayer_lock.try_update_object(player_id, new_object)
    assert success


def test_update_object_locked(test_object):
    multiplayer_lock = MultiplayerObjectLock(test_object)
    player_id = "1"
    multiplayer_lock.try_lock(player_id)
    new_object = object()
    success = multiplayer_lock.try_update_object("2", new_object)
    assert not success
