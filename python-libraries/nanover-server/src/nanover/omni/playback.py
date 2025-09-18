from os import PathLike
from pathlib import Path
from typing import List, Tuple, Iterable, Set

from nanover.app.types import AppServer
from nanover.trajectory import FrameData
from nanover.utilities.change_buffers import DictionaryChange
from nanover.recording.reading import iter_recording_files

MICROSECONDS_TO_SECONDS = 1 / 1000000
SCENE_POSE_IDENTITY = [0, 0, 0, 0, 0, 0, 1, 1, 1, 1]

Entry = Tuple[float, FrameData | None, DictionaryChange | None]


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
        traj: PathLike[str] | None = None,
        state: PathLike[str] | None = None
    ):
        self.name = name
        self.traj_path = traj
        self.state_path = state

        self.app_server: AppServer | None = None

        self.entries: List[Entry] = []
        self.changed_keys: Set[str] = set()
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
        self.changed_keys = {
            key
            for _, _, update in self.entries
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
        time, frame, update = self.entries[self.next_entry_index]
        self.next_entry_index = self.next_entry_index + 1
        self.time = time
        self.emit(frame=frame, update=update)

    def emit(self, *, frame: FrameData | None, update: DictionaryChange | None):
        if self.app_server is None:
            return

        if frame is not None:
            index = 0 if "index" not in frame.values else int(frame.values["index"])
            self.app_server.frame_publisher.send_frame(index, frame)
        if update is not None:
            self.app_server.clear_locks()
            self.app_server.update_state(None, update)
