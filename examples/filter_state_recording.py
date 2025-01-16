"""
Example of filtering recording to remove certain state keys.
"""

from typing import Set

from nanover.recording.reading import iter_recording_entries
from nanover.recording.writing import record_entries
from nanover.protocol.state import StateUpdate
from nanover.state.state_service import (
    state_update_to_dictionary_change,
    dictionary_change_to_state_update,
)
from nanover.utilities.change_buffers import DictionaryChange

in_path = "test-17ala-avatar-interacts.state"
out_path = in_path + "-filtered.state"


def is_key_unwanted(key):
    prefixes = {"imd"}
    return any(key.startswith(prefix) for prefix in prefixes)


def keys_in_change(change: DictionaryChange):
    yield from change.updates.keys()
    yield from change.removals


def filter_state_update(message: StateUpdate, keys: Set):
    change = state_update_to_dictionary_change(message)
    change.updates = {
        key: value for key, value in change.updates.items() if key not in keys
    }
    change.removals = {key for key in change.removals if key not in keys}
    message = dictionary_change_to_state_update(change)
    return message


with open(in_path, "rb") as infile:
    state_entries = list(iter_recording_entries(infile, StateUpdate))

changed_keys = {
    key
    for timestamp, change in state_entries
    for key in keys_in_change(state_update_to_dictionary_change(change))
}

print("inspecting", in_path, "\n")
print("STATE KEYS:\n" + "\n".join(sorted(changed_keys)) + "\n")

keys_to_remove = {key for key in changed_keys if is_key_unwanted(key)}

print("TO REMOVE:\n" + "\n".join(sorted(keys_to_remove)) + "\n")
print("writing", out_path, "\n")

filtered = [
    (timestamp, filter_state_update(message, keys=keys_to_remove))
    for (timestamp, message) in state_entries
]

with open(out_path, "wb") as outfile:
    record_entries(outfile, filtered)
