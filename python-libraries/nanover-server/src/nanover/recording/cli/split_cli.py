"""
Split a NanoVer recording into multiple parts.

Example when used as a cli:

.. code:: bash

    # show the help
    nanover-split-recording --help

    # split a recording whenever the frame resets
    nanover-split-recording recording.nanover.zip

    # split a recording whenever the simulation counter key changes and include a custom key's value in the resulting filenames
    nanover-split-recording recording.nanover.zip -s system.simulation.counter -n puppeteer.simulation-name

"""

import argparse
from os import PathLike
from pathlib import Path
from typing import Callable

from nanover.recording.reading import (
    RecordingEvent,
    NanoverRecordingReader,
    MessageEvent,
)
from nanover.recording.writing import NanoverRecordingWriter
from nanover.trajectory.keys import SIMULATION_COUNTER, FRAME_INDEX


def split_on_frame_reset(event: RecordingEvent):
    return (
        event.next_frame_event is not None
        and event.next_frame_event.message[FRAME_INDEX] == 0
    )


def make_key_change_predicate(key: str):
    def predicate(event: RecordingEvent):
        if event.next_frame_event is not None:
            return key in event.next_frame_event.message
        if event.next_state_event is not None:
            return (
                key in event.next_state_event.message["updates"]
                or key in event.next_state_event.message["removals"]
            )
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
    path: PathLike[str],
    *,
    split_predicate: Callable[[RecordingEvent], bool] = split_on_frame_reset,
    name_template=name_basic,
):
    """
    :param path: Path of a NanoVer recording
    :param split_predicate: Predicate function that takes information about the next frame and returns True if the
    recording should split at this point.
    :param name_template: Function that generates a filename for each recording.
    """
    input_path = Path(path)

    split_count = 0
    last_event = None
    current_base_timestamp = 0
    current_writer: NanoverRecordingWriter | None = None

    temp_path = f"{input_path.parent}/{input_path.stem}--TEMP.nanover.zip"

    def close_all():
        split_stem = name_template(
            input_stem=input_path.stem,
            index=split_count,
            last_event=last_event,
        )

        if current_writer is not None:
            current_writer.close()
            Path(temp_path).rename(input_path.parent / f"{split_stem}.nanover.zip")

    def open_all():
        nonlocal current_writer
        current_writer = NanoverRecordingWriter.from_path(temp_path)

    try:
        open_all()

        with NanoverRecordingReader.from_path(path) as reader:
            for event in reader.iter_max():
                last_event = event

                if split_predicate(event):
                    close_all()
                    split_count += 1
                    current_base_timestamp = event.timestamp
                    open_all()

                if current_writer is not None:
                    message = {}
                    if event.next_frame_event is not None:
                        message["frame"] = event.next_frame_event.message
                    if event.next_state_event is not None:
                        message["state"] = event.next_state_event.message

                    current_writer.write_message_event(
                        MessageEvent(
                            timestamp=event.timestamp - current_base_timestamp,
                            message=message,
                        )
                    )
    finally:
        close_all()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        dest="path",
        metavar="PATH",
        help="Recording file to split",
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

    if args.split_by_key is not None:
        predicate = make_key_change_predicate(args.split_by_key)
    else:
        predicate = split_on_frame_reset

    if args.name_by_key is not None:
        template = make_key_name_template(args.name_by_key)
    else:
        template = name_basic

    split_recording(
        path=args.path,
        split_predicate=predicate,
        name_template=template,
    )


if __name__ == "__main__":
    main()
