"""
Split a NanoVer recording into multiple parts.

Example when used as a cli:

.. code:: bash

    # show the help
    nanover-split-recording --help

    # split a recording whenever the frame resets
    nanover-split-recording recording.traj recording.state

    # split a recording whenever the simulation counter key changes and include a custom key's value in the resulting filenames
    nanover-split-recording recording.traj recording.state -s system.simulation.counter -n puppeteer.simulation-name

"""

import argparse
from os import PathLike
from pathlib import Path
from typing import Optional, Callable, BinaryIO

from nanover.recording.utilities import RecordingEvent, iter_recording_max
from nanover.recording.writing import write_header, write_entry
from nanover.state.state_service import dictionary_change_to_state_update
from nanover.trajectory.frame_data import SIMULATION_COUNTER
from nanover.utilities.change_buffers import DictionaryChange

from nanover.protocol.trajectory import GetFrameResponse
from nanover.utilities.protobuf_utilities import struct_to_dict


def split_on_frame_reset(event: RecordingEvent):
    return (
        event.next_frame_event is not None
        and event.next_frame_event.message.frame_index == 0
    )


def make_key_change_predicate(key: str):
    def predicate(event: RecordingEvent):
        if event.next_frame_event is not None:
            return (
                key in event.next_frame_event.message.frame.values
                or key in event.next_frame_event.message.frame.arrays
            )
        if event.next_state_event is not None:
            return key in struct_to_dict(event.next_state_event.message.changed_keys)
        return False

    return predicate


split_on_sim_counter_change = make_key_change_predicate(SIMULATION_COUNTER)


def name_basic(
    *,
    input_stem: str,
    index: int,
    last_event: RecordingEvent,
):
    return f"{input_stem}-{index:02d}"


def make_key_name_template(key: str):
    def template(
        *,
        input_stem: str,
        index: int,
        last_event: RecordingEvent,
    ):
        name = name_basic(input_stem=input_stem, index=index, last_event=last_event)
        value = get_value(last_event, key)

        return f"{name}-{value}"

    return template


def get_value(event: RecordingEvent, key: str):
    if key in event.prev_frame.value_keys:
        return event.prev_frame.values[key]
    elif key in event.prev_frame.array_keys:
        return event.prev_frame.arrays[key]
    elif key in event.prev_state:
        return event.prev_state[key]
    return None


def split_recording(
    *,
    traj: Optional[PathLike[str]] = None,
    state: Optional[PathLike[str]] = None,
    split_predicate: Callable[[RecordingEvent], bool] = split_on_frame_reset,
    name_template=name_basic,
):
    """
    :param traj: Path of a NanoVer trajectory recording
    :param state: Path of a NanoVer state recording
    :param split_predicate: Predicate function that takes information about the next frame and returns True if the
    recording should split at this point.
    :param name_template: Function that generates a filename for each recording.
    """
    split_count = 0
    last_event = None
    current_base_timestamp = 0
    current_traj_out: Optional[BinaryIO] = None
    current_state_out: Optional[BinaryIO] = None

    input_path = Path(traj if traj is not None else state)
    temp_stem = f"{input_path.stem}--TEMP"

    def close_all():
        split_stem = name_template(
            input_stem=input_path.stem,
            index=split_count,
            last_event=last_event,
        )

        if current_traj_out is not None:
            current_traj_out.close()
            Path(f"{temp_stem}.traj").rename(input_path.parent / f"{split_stem}.traj")
        if current_state_out is not None:
            current_state_out.close()
            Path(f"{temp_stem}.state").rename(input_path.parent / f"{split_stem}.state")

    def open_all():
        nonlocal current_traj_out, current_state_out
        if traj is not None:
            current_traj_out = open(f"{temp_stem}.traj", "wb")
            write_header(current_traj_out)
        if state is not None:
            current_state_out = open(f"{temp_stem}.state", "wb")
            write_header(current_state_out)

    def dict_to_state_update(dict):
        change = DictionaryChange(updates=dict)
        state = dictionary_change_to_state_update(change)
        return state

    try:
        open_all()

        for event in iter_recording_max(traj=traj, state=state):
            last_event = event

            if split_predicate(event):
                # end the previous recording
                close_all()

                # start the next recording
                split_count += 1
                current_base_timestamp = event.timestamp
                open_all()

                # start the next recording with the outcome of the current event
                if current_traj_out is not None:
                    write_entry(
                        current_traj_out,
                        0,
                        GetFrameResponse(frame_index=0, frame=event.next_frame.raw),
                    )
                if current_state_out is not None:
                    write_entry(
                        current_state_out, 0, dict_to_state_update(event.next_state)
                    )
            else:
                timestamp = event.timestamp - current_base_timestamp
                if event.next_frame_event is not None:
                    write_entry(
                        current_traj_out, timestamp, event.next_frame_event.message
                    )
                if event.next_state_event is not None:
                    write_entry(
                        current_state_out, timestamp, event.next_state_event.message
                    )
    finally:
        close_all()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        dest="paths",
        nargs="+",
        metavar="PATH",
        help="Recording files to split (one or both of .traj and .state)",
    )
    parser.add_argument(
        "-s",
        "--split-by-key",
        help="Split when a message changes the value of (or removes) this key",
        type=str,
    )
    parser.add_argument(
        "-n",
        "--name-by-key",
        help="Include the value of this key as part of the filename of output recording",
        type=str,
    )
    args = parser.parse_args()

    paths = [Path(path) for path in args.paths]
    traj_path = next((path for path in paths if path.suffix == ".traj"), None)
    state_path = next((path for path in paths if path.suffix == ".state"), None)

    if args.split_by_key is not None:
        predicate = make_key_change_predicate(args.split_by_key)
    else:
        predicate = split_on_frame_reset

    if args.name_by_key is not None:
        template = make_key_name_template(args.name_by_key)
    else:
        template = name_basic

    split_recording(
        traj=traj_path,
        state=state_path,
        split_predicate=predicate,
        name_template=template,
    )


if __name__ == "__main__":
    main()
