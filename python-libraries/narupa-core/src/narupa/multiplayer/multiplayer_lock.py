# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.

"""
Module providing the classes for shared editing of a multiplayer object.
"""
import datetime
import threading
from threading import Lock


class MultiplayerObjectLock(object):
    """
    Used to facilitate read/write access to an object that is being shared over a multiplayer session.
    A client can request to lock an object using their player ID, during which time they alone can access the object.

    Note that the object is still read/writeable within the controlling process.


    :param max_lock_time: Maximum lock time, in seconds.

    """

    _locked: bool
    _threadlock: Lock
    _timelocked: datetime.datetime
    _max_lock_time: datetime.timedelta

    def __init__(self, locked_object, max_lock_time=20):

        self._threadlock = threading.Lock()
        self._unlock()
        self._max_lock_time = datetime.timedelta(seconds=max_lock_time)
        self._locked_object = locked_object

    @property
    def locked_object(self):
        """
        Returns the current state of the locked object.
        Note that the result should not be edited outside of the management of the owner lock,
        as that defeats the point of the lock.
        :return: The locked object.
        """
        return self._locked_object

    def is_locked(self):
        """
        Determines whether the object is currently locked.
        :return:
        """
        return self._player_id is not None and not self._lock_expired()

    def try_update_object(self, player_id, object_update):
        if self.try_lock(player_id):
            self._locked_object = object_update
            return True
        else:
            return False

    def try_lock(self, player_id) -> bool:
        """
        Try to lock the object with the given player_id.
        :param player_id: a unique identifier representing the user attempting to lock the object.
        :return: True if successfully locked, false otherwise.
        """
        with self._threadlock:
            if self.is_locked():
                return self.is_lock_owner(player_id)
            self._lock(player_id)
            return True

    def try_unlock(self, player_id):
        """
        Try to unlock the object with the given player_id.
        :param player_id: token representing lock ownership
        :return: True if successfully unlocked, false otherwise.
        """
        with self._threadlock:
            if not self.is_locked():
                return True
            if self._lock_expired() or self.is_lock_owner(player_id):
                self._unlock()
                return True
            return False

    def is_lock_owner(self, player_id) -> bool:
        """
        Determine whether the given player_id is the owner of the current lock.
        :param player_id:
        :return: If the lock is locked, returns whether the player_id matches the owner of the locker, otherwise, False.
        """
        if not self.is_locked():
            return False
        return self._player_id == player_id

    def get_lock_owner(self):
        """
        Get the locker ID.
        :return:
        """
        return self._player_id

    def _lock_expired(self) -> bool:
        now = datetime.datetime.now()
        if self._timelocked + self._max_lock_time < now:
            return True
        return False

    def _lock(self, player_id):
        self._timelocked = datetime.datetime.now()
        self._player_id = player_id

    def _unlock(self):
        self._player_id = None
        self._timelocked = None
