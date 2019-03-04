from threading import Lock
import threading
import datetime
import time

class OwnerLock(object):
    """
    Represents a lock based on ownership by guid. 
    This does not actually lock in the sense of thread locking, but can be used to control access to some property.
    """

    _locked: bool
    _guid: int
    _threadlock: Lock
    _timelocked: datetime.datetime
    _max_lock_time: datetime.timedelta

    def __init__(self, max_lock_time=20):
        """
        Initialises the owner lock.
        :param max_lock_time: Maximum lock time, in seconds.
        """
        self._threadlock = threading.Lock()
        self._unlock()
        self._max_lock_time = datetime.timedelta(seconds = max_lock_time)

    def _lock_expired(self) -> bool:
        if self._timelocked is None:
            return True
        now = datetime.datetime.now()
        if self._timelocked + self._max_lock_time < now:
            return True
        return False

    def try_lock(self, guid) -> bool:
        """
        Try to lock the scene properties with the given guid.
        :param guid: token representing lock.
        :return: True if successfully locked, false otherwise.
        """
        with self._threadlock:
            if self._locked:
                if self._lock_expired():
                    self._lock(guid)
                    return True
                else:
                    return self.is_locker(guid)
            self._lock(guid)
            return True

    def is_locked(self):
        return self._locked

    def is_locker(self, guid)-> bool:
        """
        Determine whether the given guid is the owner of the current lock.
        :param guid: 
        :return: If the lock is locked, returns whether the guid matches the owner of the locker, otherwise, False. 
        """
        if self._locked == False: return False
        return self._guid == guid

    def _lock(self, guid):
        self._timelocked = datetime.datetime.now()
        self._locked = True
        self._guid = guid

    def _unlock(self):
        self._locked = False
        self._guid = 0
        self._timelocked = None

    def try_unlock(self, guid):
        """
        Try to unlock the scene properties with the given guid.
        :param guid: token representing lock ownership
        :return: True if successfully unlocked, false otherwise.
        """
        with self._threadlock:
            if self._locked is False:
                return True
            if self._lock_expired() or self.is_locker(guid):
                self._unlock()
                return True
            return False
