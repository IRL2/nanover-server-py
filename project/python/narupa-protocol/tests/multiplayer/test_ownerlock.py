from time import sleep

from narupy.multiplayer.scene import OwnerLock

def test_lock():
    owner_lock = OwnerLock()
    guid = 1
    result = owner_lock.try_lock(guid)
    assert result == True
    assert owner_lock.is_locked() == True

def test_lock_ownership():
    owner_lock = OwnerLock()
    guid = 1
    owner_lock.try_lock(guid)
    result = owner_lock.is_locker(guid)
    assert result == True

def test_unlock_unlocked():
    owner_lock = OwnerLock()
    guid = 1
    result = owner_lock.try_unlock(guid)
    assert result == True

def test_unlock():
    owner_lock = OwnerLock()
    guid = 1
    owner_lock.try_lock(guid)
    result = owner_lock.try_unlock(guid)
    assert result == True
    assert owner_lock.is_locked() == False

def test_second_lock_attempt():
    owner_lock = OwnerLock()
    guid = 1
    owner_lock.try_lock(guid)

    guid = 2
    result = owner_lock.try_lock(guid)
    assert result == False

def test_second_unlock_attempt():
    owner_lock = OwnerLock()
    guid = 1
    owner_lock.try_lock(guid)

    guid = 2
    result = owner_lock.try_unlock(guid)
    assert result == False

def test_elapsed_lock():
    owner_lock = OwnerLock(max_lock_time=1)
    guid = 1
    owner_lock.try_lock(guid)

    sleep(1)

    guid = 2
    result = owner_lock.try_lock(guid)
    assert result == True

def test_elapsed_unlock():
    owner_lock = OwnerLock(max_lock_time=1)
    guid = 1
    owner_lock.try_lock(guid)

    sleep(1)

    guid = 2
    result = owner_lock.try_unlock(guid)
    assert result == True

