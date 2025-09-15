"""
Check a NanoVer recording for errors.

Example when used as a cli:

.. code:: bash

    # show the help
    nanover-check-recording --help

    # check a recording
    nanover-check-recording recording.traj recording.state

"""

import argparse
from contextlib import suppress
from dataclasses import dataclass
from os import PathLike, SEEK_END, SEEK_CUR
from pathlib import Path
from typing import Optional

from nanover.recording.reading import (
    read_header,
    MessageRecordingReader,
    read_entry_info,
)
from nanover.recording.writing import write_header, write_buffer


@dataclass(kw_only=True)
class IndexRepairEntry:
    offset: int
    timestamp: float
    size: int


@dataclass(kw_only=True)
class IndexRepairDefectiveBlock:
    offset: int
    size: int
    prev: Optional[IndexRepairEntry] = None
    next: Optional[IndexRepairEntry] = None


def reindex_repair(
    reader: MessageRecordingReader,
    max_interval=120,
):
    reader.io.seek(0, SEEK_END)
    max_seek = reader.io.tell()
    max_interval_ms = max_interval * 10**6

    entries: list[IndexRepairEntry] = []
    defective_blocks: list[IndexRepairDefectiveBlock] = []

    def add_entry(entry: IndexRepairEntry):
        entries.append(entry)

    def is_valid(entry: IndexRepairEntry):
        if not entries:
            return True

        prev = entries[-1]

        monotonic = entry.timestamp > prev.timestamp
        continuous = entry.timestamp < prev.timestamp + max_interval_ms
        bounded = entry.offset + entry.size < max_seek

        return monotonic and continuous and bounded

    reader.io.seek(0)
    read_header(reader.io)

    with suppress(EOFError):
        defective = None

        while True:
            offset = reader.io.tell()
            timestamp, size = read_entry_info(reader.io)

            entry = IndexRepairEntry(
                offset=offset,
                timestamp=timestamp,
                size=size,
            )

            if is_valid(entry):
                if defective is not None:
                    defective.next = entry
                add_entry(entry)
                defective = None
                reader.io.seek(size, SEEK_CUR)
            else:
                if defective is None:
                    defective = IndexRepairDefectiveBlock(
                        offset=offset, size=0, prev=entries[-1]
                    )
                    defective_blocks.append(defective)
                defective.size += 1
                reader.io.seek(1, SEEK_CUR)

    reader.message_offsets = [entry.offset for entry in entries]
    return entries, defective_blocks


def check_recording(path: PathLike[str]):
    path = Path(path)

    with open(path, "rb") as infile:
        reader = MessageRecordingReader(infile)
        entries, defective_blocks = reindex_repair(reader)

        for defective in defective_blocks:
            prev_time = (
                defective.prev.timestamp / 10**6
                if defective.prev is not None
                else "START"
            )
            next_time = (
                defective.next.timestamp / 10**6
                if defective.next is not None
                else "END"
            )
            print(
                f"DEFECTIVE BLOCK ({defective.size}) between {prev_time:.2f}s and {next_time:.2f}s"
            )

        if defective_blocks:
            outpath = f"{path.with_suffix("")}--repaired--{path.suffix}"

            with open(outpath, "wb") as outfile:
                write_header(outfile)
                for entry in reader:
                    write_buffer(outfile, entry.timestamp, entry.buffer)

            print(f"Wrote repaired file to {outpath}")
        else:
            print(f"No defects in {path}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        dest="paths",
        nargs="+",
        metavar="PATH",
        help="Recording files to check (one or both of .traj and .state)",
    )
    args = parser.parse_args()

    paths = [Path(path) for path in args.paths]
    traj_path = next((path for path in paths if path.suffix == ".traj"), None)
    state_path = next((path for path in paths if path.suffix == ".state"), None)

    if traj_path is not None:
        check_recording(traj_path)
    if state_path is not None:
        check_recording(state_path)


if __name__ == "__main__":
    main()
