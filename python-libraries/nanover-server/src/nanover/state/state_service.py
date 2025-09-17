"""
Module providing an implementation of the :class:`StateServicer`.
"""

from typing import Set, Dict, ContextManager
from nanover.utilities.protobuf_utilities import (
    deep_copy_serializable_dict,
    struct_to_dict,
    dict_to_struct,
    Serializable,
)
from nanover.utilities.change_buffers import (
    DictionaryChange,
    DictionaryChangeBuffer,
)
from nanover.protocol.state import (
    StateUpdate,
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

    def lock_state(self) -> ContextManager[Dict[str, Serializable]]:
        """
        Context manager for reading the current state while delaying any changes
        to it.
        """
        return self.state_dictionary.lock_content()

    def copy_state(self) -> Dict[str, Serializable]:
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
        :raises TypeError: if the update values cannot be serialized for
            transmission.
        """
        validate_dict_is_serializable(change.updates)
        self.state_dictionary.update_state(access_token, change)

    def update_locks(
        self,
        access_token: Serializable,
        acquire: Dict[str, float] | None = None,
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


def validate_dict_is_serializable(dictionary):
    """
    :raises TypeError: if the given dictionary cannot be converted to a protobuf
        struct.
    """
    try:
        dict_to_struct(dictionary)
    except ValueError as e:
        raise TypeError("Data is not serializable with protobuf.") from e


def state_update_to_dictionary_change(update: StateUpdate) -> DictionaryChange:
    """
    Convert a protobuf StateUpdate to a DictionaryChange.

    :param update: a protobuf StateUpdate which encodes key removals as keys
        with a protobuf null value.
    :return: an equivalent DictionaryChange representing the key changes and
        key removals of the StateUpdate.
    """
    removals = set()
    updates = {}

    for key, value in struct_to_dict(update.changed_keys).items():
        if value is None:
            removals.add(key)
        else:
            updates[key] = value

    return DictionaryChange(removals=removals, updates=updates)


def dictionary_change_to_state_update(change: DictionaryChange) -> StateUpdate:
    """
    Convert a DictionaryChange to a protobuf StateUpdate.

    :param change: a DictionaryChange which species key changes and key removals
        to make to a dictionary.
    :return: an equivalent protobuf StateUpdate representing the key removals
        as key changes to a protobuf null value.
    """
    changes, removals = change

    update = StateUpdate()
    update.changed_keys.update(changes)
    for key in removals:
        update.changed_keys[key] = None

    return update
