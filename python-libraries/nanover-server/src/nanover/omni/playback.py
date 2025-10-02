from os import PathLike
from pathlib import Path

from nanover.app.types import AppServer
from nanover.recording.reading import (
    RecordingIndexEntry,
    NanoverRecordingReader,
)
from nanover.trajectory import FrameData
from nanover.trajectory.keys import FRAME_INDEX
from nanover.utilities.change_buffers import DictionaryChange

MICROSECONDS_TO_SECONDS = 1 / 1000000
SCENE_POSE_IDENTITY = [0, 0, 0, 0, 0, 0, 1, 1, 1, 1]

Entry = tuple[float, RecordingIndexEntry]


class PlaybackSimulation:
    @classmethod
    def from_path(
        cls,
        *,
        name: str | None = None,
        path: PathLike[str],
    ):
        path = Path(path)
        name = name or path.stem

        return cls(name=name, path=path)

    def __init__(
        self,
        *,
        name: str,
        path: PathLike[str],
    ):
        self.name = name
        self.path = path

        self.app_server: AppServer | None = None

        self.entries: list[Entry] = []
        self.changed_keys: set[str] = set()
        self.next_entry_index = 0
        self.time = 0.0

        self.reader = NanoverRecordingReader.from_path(path)

    def load(self):
        """
        Load and set up the simulation if it isn't done already.
        """
        self.entries = [
            (entry.metadata["timestamp"] * MICROSECONDS_TO_SECONDS, entry)
            for entry in self.reader
        ]

        updates = (self.reader.get_state_from_entry(entry) for entry in self.reader)

        self.changed_keys = {
            key
            for update in updates
            if update is not None
            for key in update.updates.keys()
            if key != "scene"
        }

    def reset(self, app_server: AppServer):
        """
        Reset the playback to its initial state.
        :param app_server: The app server hosting the frame publisher and imd state
        """
        self.app_server = app_server
        self.next_entry_index = 0
        self.time = 0.0

        # clear simulation and reset box pose to identity
        self.emit(
            frame=FrameData(),
            update=DictionaryChange(
                updates={"scene": SCENE_POSE_IDENTITY}, removals=self.changed_keys
            ),
        )

    def advance_by_one_step(self):
        """
        Advance playback to the next point a frame or update should be reported, and report it.
        """
        assert self.app_server is not None
        try:
            self.advance_to_next_entry()
        except IndexError:
            self.reset(self.app_server)
            self.advance_to_next_entry()

    def advance_by_seconds(self, dt: float):
        """
        Advance the playback by some seconds, emitting any intermediate frames and state updates.
        :param dt: Time to advance playback by in seconds
        """
        assert self.app_server is not None
        next_time = self.time + dt

        try:
            while self.entries[self.next_entry_index][0] <= next_time:
                self.advance_to_next_entry()
        except IndexError:
            self.reset(self.app_server)
        else:
            self.time = next_time

    def advance_to_next_entry(self):
        """
        Advance playback to the next point a frame or update should be reported, and report it.
        """
        self.time, entry = self.entries[self.next_entry_index]
        self.next_entry_index = self.next_entry_index + 1

        self.emit(
            frame=self.reader.get_frame_from_entry(entry),
            update=self.reader.get_state_from_entry(entry),
        )

    def emit(self, *, frame: FrameData | None, update: DictionaryChange | None):
        if self.app_server is None:
            return

        if frame is not None:
            index = 0 if FRAME_INDEX not in frame else frame.frame_index
            self.app_server.frame_publisher.send_frame(index, frame)
        if update is not None:
            self.app_server.clear_locks()
            self.app_server.update_state(None, update)
