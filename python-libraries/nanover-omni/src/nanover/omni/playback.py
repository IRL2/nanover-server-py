from contextlib import suppress
from os import PathLike
from pathlib import Path
from typing import List, Optional, Tuple, Iterable

from nanover.app import NanoverImdApplication
from nanover.trajectory import FrameData
from nanover.utilities.change_buffers import DictionaryChange
from nanover.recording.reading import iter_recording_files
from nanover.utilities.key_lockable_map import ResourceLockedError

MICROSECONDS_TO_SECONDS = 1 / 1000000

Entry = Tuple[float, Optional[FrameData], Optional[DictionaryChange]]


class PlaybackSimulation:
    @classmethod
    def from_paths(cls, paths: Iterable[PathLike[str]]):
        """
        Construct this from one or both of trajectory and state recording file paths.
        :param paths: One or both of trajectory and state recording file paths.
        """
        paths = [Path(path) for path in paths]
        traj_path = next((path for path in paths if path.suffix == ".traj"), None)
        state_path = next((path for path in paths if path.suffix == ".state"), None)

        return cls(
            name=paths[0].stem,
            traj=traj_path,
            state=state_path,
        )

    def __init__(
        self,
        name,
        *,
        traj: Optional[PathLike[str]] = None,
        state: Optional[PathLike[str]] = None
    ):
        self.name = name
        self.traj_path = traj
        self.state_path = state

        self.app_server: Optional[NanoverImdApplication] = None

        self.entries: List[Entry] = []
        self.next_entry_index = 0
        self.time = 0.0

    def load(self):
        """
        Load and set up the simulation if it isn't done already.
        """
        entries = iter_recording_files(traj=self.traj_path, state=self.state_path)
        self.entries = [
            (time * MICROSECONDS_TO_SECONDS, frame, update)
            for time, frame, update in entries
        ]

    def reset(self, app_server: NanoverImdApplication):
        """
        Reset the playback to its initial state.
        :param app_server: The app server hosting the frame publisher and imd state
        """
        self.app_server = app_server
        self.next_entry_index = 0
        self.time = 0.0

        # clear simulation
        self.emit(frame=FrameData(), update=None)

    def advance_by_one_step(self):
        """
        Advance playback to the next point a frame or update should be reported, and report it.
        """
        try:
            self.advance_to_next_entry()
        except IndexError:
            self.next_entry_index = 0
            self.time = 0.0
            self.advance_to_next_entry()

    def advance_by_seconds(self, dt: float):
        """
        Advance the playback by some seconds, emitting any intermediate frames and state updates.
        :param dt: Time to advance playback by in seconds
        """
        next_time = self.time + dt

        try:
            while self.entries[self.next_entry_index][0] <= next_time:
                self.advance_to_next_entry()
        except IndexError:
            self.next_entry_index = 0
            self.time = 0.0
        else:
            self.time = next_time

    def advance_to_next_entry(self):
        """
        Advance playback to the next point a frame or update should be reported, and report it.
        """
        time, frame, update = self.entries[self.next_entry_index]
        self.next_entry_index = self.next_entry_index + 1
        self.time = time
        self.emit(frame=frame, update=update)

    def emit(self, *, frame: Optional[FrameData], update: Optional[DictionaryChange]):
        if self.app_server is None:
            return

        if frame is not None:
            index = 0 if "index" not in frame.values else int(frame.values["index"])
            self.app_server.frame_publisher.send_frame(index, frame)
        if update is not None:
            with suppress(ResourceLockedError):
                self.app_server.server.update_state(None, update)
