from dataclasses import dataclass, field
from itertools import groupby
from os import PathLike
from typing import Any, IO, BinaryIO
from zipfile import ZipFile

import msgpack

from nanover.utilities.state_dictionary import StateDictionary
from nanover.trajectory import FrameData
from nanover.trajectory.convert import (
    unpack_dict_frame,
    convert_dict_state_to_dictionary_change,
)
from nanover.trajectory.keys import SIMULATION_COUNTER
from nanover.utilities.change_buffers import DictionaryChange

RECORDING_INDEX_FILENAME = "index.msgpack"
RECORDING_MESSAGES_FILENAME = "messages.msgpack"


@dataclass(kw_only=True)
class MessageEvent:
    timestamp: int | None
    message: dict[str, Any]


@dataclass(kw_only=True)
class FrameRecordingEvent(MessageEvent):
    prev_frame: FrameData
    next_frame: FrameData

    @classmethod
    def make_empty(cls):
        return cls(
            timestamp=0,
            message={},
            prev_frame=FrameData(),
            next_frame=FrameData(),
        )


@dataclass(kw_only=True)
class StateRecordingEvent(MessageEvent):
    prev_state: dict
    next_state: dict

    @classmethod
    def make_empty(cls):
        return cls(
            timestamp=0,
            message={},
            prev_state={},
            next_state={},
        )


@dataclass(kw_only=True)
class RecordingEvent:
    prev_frame_event: FrameRecordingEvent
    prev_state_event: StateRecordingEvent
    next_frame_event: FrameRecordingEvent | None
    next_state_event: StateRecordingEvent | None

    @property
    def timestamp(self) -> int | None:
        if self.next_frame_event is not None:
            return self.next_frame_event.timestamp
        elif self.next_state_event is not None:
            return self.next_state_event.timestamp
        return None

    @property
    def prev_frame(self):
        return self.prev_frame_event.next_frame

    @property
    def next_frame(self):
        return (
            self.next_frame_event.next_frame
            if self.next_frame_event is not None
            else self.prev_frame
        )

    @property
    def prev_state(self):
        return self.prev_state_event.next_state

    @property
    def next_state(self):
        return (
            self.next_state_event.next_state
            if self.next_state_event is not None
            else self.prev_state
        )


@dataclass(kw_only=True)
class RecordingIndexEntry:
    offset: int
    length: int
    metadata: dict[str, Any] = field(default_factory=dict)

    def read_from(self, io: IO[bytes]):
        io.seek(self.offset)
        return io.read(self.length)


class MessageZipReader:
    @classmethod
    def from_path(cls, path: PathLike[str]):
        """
        Read a recording from a filepath.
        """
        return cls(ZipFile(path, "r"))

    @classmethod
    def from_io(cls, io: BinaryIO):
        """
        Read a recording from binary data source.
        """
        return cls(ZipFile(io, "r"))

    def __init__(self, zipfile: ZipFile):
        self.zipfile = zipfile
        self.messagesfile = zipfile.open(RECORDING_MESSAGES_FILENAME)
        self.index = parse_index(zipfile)

    def close(self):
        self.zipfile.close()

    def get_message_from_entry(self, entry: RecordingIndexEntry) -> dict[str, Any]:
        return msgpack.unpackb(entry.read_from(self.messagesfile))

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def __iter__(self):
        yield from self.index

    def __len__(self):
        return len(self.index)

    def __getitem__(self, index):
        return self.index[index]


class NanoverRecordingReader(MessageZipReader):
    def iter_max(self):
        """
        Iterate recording yielding recording events in timestamp order, with each event containing the full
        information of previous frame, previous state, current message, next frame, and next state.
        """
        current_frame = FrameData()
        current_state = StateDictionary()

        prev_frame_event = FrameRecordingEvent.make_empty()
        prev_state_event = StateRecordingEvent.make_empty()

        for entry in self.index:
            next_frame_event: FrameRecordingEvent | None = None
            next_state_event: StateRecordingEvent | None = None

            timestamp: int | None = entry.metadata.get("timestamp", None)
            frame_message = self.get_message_type_from_entry(entry, "frame")
            state_message = self.get_message_type_from_entry(entry, "state")

            if frame_message is not None:
                prev_frame = current_frame.copy()
                current_frame.update(FrameData(frame_message))

                next_frame_event = FrameRecordingEvent(
                    timestamp=timestamp,
                    message=frame_message,
                    prev_frame=prev_frame,
                    next_frame=current_frame,
                )
            if state_message is not None:
                prev_state = current_state.copy_content()
                current_state.update_state(
                    None, convert_dict_state_to_dictionary_change(state_message)
                )

                next_state_event = StateRecordingEvent(
                    timestamp=timestamp,
                    message=state_message,
                    prev_state=prev_state,
                    next_state=current_state.copy_content(),
                )

            yield RecordingEvent(
                prev_frame_event=prev_frame_event,
                prev_state_event=prev_state_event,
                next_frame_event=next_frame_event,
                next_state_event=next_state_event,
            )

            prev_frame_event = next_frame_event or prev_frame_event
            prev_state_event = next_state_event or prev_state_event

    def iter_frames_full(self):
        current = FrameData()
        for entry, frame in self.iter_frame_updates():
            current.update(frame)
            yield entry, current.copy()

    def iter_states_full(self):
        current = StateDictionary()
        for entry, change in self.iter_state_updates():
            current.update_state(None, change)
            yield entry, current.copy_content()

    def iter_frame_updates(self):
        for entry in self:
            frame = self.get_frame_from_entry(entry)
            if frame is not None:
                yield entry, frame

    def iter_state_updates(self):
        for entry in self:
            state = self.get_state_from_entry(entry)
            if state is not None:
                yield entry, state

    def iter_message_events(self):
        for entry in self:
            message = self.get_message_from_entry(entry)
            yield MessageEvent(
                timestamp=entry.metadata["timestamp"],
                message=message,
            )

    def get_frame_from_entry(self, entry: RecordingIndexEntry) -> FrameData | None:
        if "frame" not in entry.metadata.get("types", ("frame",)):
            return None

        message = self.get_message_from_entry(entry)

        if "frame" not in message:
            return None

        return FrameData(unpack_dict_frame(message["frame"]))

    def get_state_from_entry(
        self, entry: RecordingIndexEntry
    ) -> DictionaryChange | None:
        if "state" not in entry.metadata.get("types", ("state",)):
            return None

        message = self.get_message_from_entry(entry)

        if "state" not in message:
            return None

        return convert_dict_state_to_dictionary_change(message["state"])

    def get_message_type_from_entry(
        self, entry: RecordingIndexEntry, type: str
    ) -> Any | None:
        if type not in entry.metadata.get("types", (type,)):
            return None
        message = self.get_message_from_entry(entry)
        return message.get(type, None)


def split_by_simulation_counter(path: PathLike[str]):
    def get_simulation_counter(triplet):
        _, frame, _ = triplet
        return frame.frame_dict.get(SIMULATION_COUNTER)

    full_view = iter_full_view(path)
    sessions = [
        list(session)
        for simulation_counter, session in groupby(
            full_view, key=get_simulation_counter
        )
        if simulation_counter is not None
    ]

    return sessions


def iter_full_view(path: PathLike[str]):
    full_frame = FrameData()
    full_state = StateDictionary()

    for time, frame, update in iter_recording_file(path):
        if frame is not None:
            full_frame.update(frame)

        if update is not None:
            full_state.update_state(None, update)

        yield time, full_frame.copy(), full_state.copy_content()


def iter_recording_file(path: PathLike[str]):
    with MessageZipReader.from_path(path) as reader:
        for entry in reader:
            frame, update = None, None
            message = reader.get_message_from_entry(entry)
            if "frame" in message:
                frame = FrameData(unpack_dict_frame(message["frame"]))
            if "state" in message:
                update = convert_dict_state_to_dictionary_change(message["state"])

            if frame is not None or update is not None:
                yield entry.metadata["timestamp"], frame, update


def message_events_from_recording(path: PathLike[str]):
    with MessageZipReader.from_path(path) as reader:
        for entry in reader:
            yield MessageEvent(
                timestamp=entry.metadata["timestamp"],
                message=reader.get_message_from_entry(entry),
            )


def parse_index(zipfile: ZipFile):
    with zipfile.open(RECORDING_INDEX_FILENAME) as index_file:
        return [
            RecordingIndexEntry(**entry) for entry in msgpack.unpackb(index_file.read())
        ]
