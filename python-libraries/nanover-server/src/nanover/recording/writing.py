from dataclasses import asdict
from io import BytesIO
from os import PathLike
from typing import Iterable, Any, BinaryIO
from zipfile import ZipFile

import msgpack
import numpy as np

from nanover.recording.reading import (
    RECORDING_INDEX_FILENAME,
    RECORDING_MESSAGES_FILENAME,
    MessageEvent,
    RecordingIndexEntry,
)
from nanover.trajectory.keys import FRAME_INDEX


class NanoverRecordingWriter:
    @classmethod
    def from_path(cls, path: PathLike[str] | str):
        return cls(path)

    def __init__(self, outfile: BinaryIO | PathLike[str] | str):
        self.offset = 0
        self.archive = ZipFile(outfile, "w")
        self.messages_file = self.archive.open(
            RECORDING_MESSAGES_FILENAME, "w", force_zip64=True
        )
        self.index: list[RecordingIndexEntry] = []

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def close(self):
        self.messages_file.close()
        self.write_index()
        self.archive.close()

    def write_message_event(self, event: MessageEvent):
        data = msgpack.packb(event.message, default=_fallback_encoder)
        metadata = generate_metadata(event.message)
        metadata["timestamp"] = event.timestamp
        entry = RecordingIndexEntry(
            offset=self.offset,
            length=len(data),
            metadata=metadata,
        )
        self.index.append(entry)
        self.messages_file.write(data)
        self.offset += entry.length

    def write_index(self):
        with self.archive.open(RECORDING_INDEX_FILENAME, "w") as index_file:
            data = [asdict(entry) for entry in self.index]
            index_file.write(msgpack.packb(data))


def record_messages(outfile: str | BytesIO, messages: Iterable[MessageEvent]):
    with NanoverRecordingWriter(outfile) as writer:
        for event in messages:
            writer.write_message_event(event)


def generate_metadata(message: dict[str, Any]) -> dict[str, Any]:
    metadata = {
        "types": list(message.keys()),
    }

    if "frame" in message and FRAME_INDEX in message["frame"]:
        metadata[FRAME_INDEX] = message["frame"][FRAME_INDEX]

    return metadata


def _fallback_encoder(obj: Any) -> Any:
    """
    Converts, if possible, a type msgpack doesn't understand into a basic type it can encode.

    :param obj: object to be converted
    :return: simplified object
    """
    # encode numpy arrays as simple lists
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    raise TypeError(f"Unknown type: {obj}")
