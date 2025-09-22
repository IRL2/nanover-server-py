from dataclasses import asdict
from typing import Iterable, Any
from zipfile import ZipFile

import msgpack

from nanover.recording2.reading import (
    RECORDING_INDEX_FILENAME,
    RECORDING_MESSAGES_FILENAME,
    MessageEvent,
    RecordingIndexEntry,
)
from nanover.trajectory.frame_data import FRAME_INDEX


def record_messages(out_path: str, messages: Iterable[MessageEvent]):
    with ZipFile(out_path, "w") as archive:
        entries: list[RecordingIndexEntry] = []
        try:
            with archive.open(
                RECORDING_MESSAGES_FILENAME, "w", force_zip64=True
            ) as messages_file:
                offset = 0
                for event in messages:
                    data = msgpack.packb(event.message)
                    metadata = generate_metadata(event.message)
                    metadata["timestamp"] = event.timestamp
                    entry = RecordingIndexEntry(
                        offset=offset,
                        length=len(data),
                        metadata=metadata,
                    )
                    entries.append(entry)
                    messages_file.write(data)
                    offset += entry.length
        finally:
            with archive.open(RECORDING_INDEX_FILENAME, "w") as index_file:
                data = [asdict(entry) for entry in entries]
                index_file.write(msgpack.packb(data))


def generate_metadata(message: dict[str, Any]) -> dict[str, Any]:
    metadata = {
        "types": list(message.keys()),
    }

    if "frame" in message and FRAME_INDEX in message["frame"]:
        metadata[FRAME_INDEX] = message["frame"][FRAME_INDEX]

    return metadata
