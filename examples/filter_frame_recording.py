"""
Example of filtering recording to remove certain frame keys.
"""

from typing import Set

from nanover.recording.reading import iter_recording_entries
from nanover.recording.writing import record_entries
from nanover.protocol.trajectory import GetFrameResponse
from nanover.trajectory import FrameData

in_path = "test-17ala-avatar-interacts.traj"
out_path = in_path + "-filtered.traj"


def is_key_unwanted(key):
    prefixes = {"system.simulation.counter"}
    return any(key.startswith(prefix) for prefix in prefixes)


def keys_in_frame(frame_data: FrameData):
    yield from frame_data.array_keys
    yield from frame_data.value_keys


def filter_frame(message: GetFrameResponse, keys: Set):
    frame = FrameData(message.frame)

    for key in keys:
        del frame[key]

    return message


with open(in_path, "rb") as infile:
    frame_entries = list(iter_recording_entries(infile, GetFrameResponse))

changed_keys = {
    key
    for timestamp, message in frame_entries
    for key in keys_in_frame(FrameData(message.frame))
}

print("inspecting", in_path, "\n")
print("FRAME KEYS:\n" + "\n".join(sorted(changed_keys)) + "\n")

keys_to_remove = {key for key in changed_keys if is_key_unwanted(key)}

print("TO REMOVE:\n" + "\n".join(sorted(keys_to_remove)) + "\n")
print("writing", out_path, "\n")

filtered = [
    (timestamp, filter_frame(message, keys=keys_to_remove))
    for (timestamp, message) in frame_entries
]

with open(out_path, "wb") as outfile:
    record_entries(outfile, filtered)
