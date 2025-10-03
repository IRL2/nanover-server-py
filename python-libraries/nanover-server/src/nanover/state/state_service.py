"""
Module providing an implementation of the :class:`StateServicer`.
"""

from typing import Set, ContextManager
from nanover.utilities.protobuf_utilities import (
    deep_copy_serializable_dict,
    Serializable,
)
from nanover.utilities.change_buffers import (
    DictionaryChange,
    DictionaryChangeBuffer,
)
from .state_dictionary import StateDictionary


class StateService:
    """
    Implementation of the State service, for tracking and making changes to a
    shared key/value store.
    """

    state_dictionary: StateDictionary

    def __init__(self):
        super().__init__()
        self.state_dictionary = StateDictionary()
        self.name: str = "service"
        self._id = "service"

    def close(self):
        self.state_dictionary.freeze()

    def lock_state(self) -> ContextManager[dict[str, Serializable]]:
        """
        Context manager for reading the current state while delaying any changes
        to it.
        """
        return self.state_dictionary.lock_content()

    def copy_state(self) -> dict[str, Serializable]:
        """
        Return a deep copy of the current state.
        """
        with self.lock_state() as state:
            return deep_copy_serializable_dict(state)

    def update_state(self, access_token: Serializable, change: DictionaryChange):
        """
        Attempts an atomic update of the shared key/value store. If any key
        cannot be updated, no change will be made.

        :raises ResourceLockedError: if the access token cannot acquire all keys
            for updating.
        """
        self.state_dictionary.update_state(access_token, change)

    def update_locks(
        self,
        access_token: Serializable,
        acquire: dict[str, float] | None = None,
        release: Set[str] | None = None,
    ):
        """
        Attempts to acquire and release locks on keys in the shared key/value
        store. If any of the locks cannot be acquired, none of them will be.
        Requested lock releases are carried out regardless.

        :raises ResourceLockedError: if the access token cannot acquire all
            requested keys.
        """
        self.state_dictionary.update_locks(access_token, acquire, release)

    def clear_locks(self):
        """
        Release all locks on all keys.
        """
        self.state_dictionary.clear_locks()

    def get_change_buffer(self) -> ContextManager[DictionaryChangeBuffer]:
        """
        Return a DictionaryChangeBuffer that tracks changes to this service's
        state.
        """
        return self.state_dictionary.get_change_buffer()
