from os import PathLike
from pathlib import Path
from typing import List, Optional, Tuple, Iterable

from nanover.app import NanoverImdApplication
from nanover.recording.parsing import (
    FrameEntry,
    iter_recording_files,
)
from nanover.trajectory import FrameData
from nanover.utilities.change_buffers import DictionaryChange

MICROSECONDS_TO_SECONDS = 1 / 1000000

Entry = Tuple[float, Optional[FrameEntry], Optional[DictionaryChange]]


class PlaybackSimulation:
    @classmethod
    def from_paths(cls, paths: Iterable[PathLike[str]]):
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
        self.frame_index = 0
        self.time = 0

        self.paused = False

    def load(self):
        entries = iter_recording_files(traj=self.traj_path, state=self.state_path)
        self.entries = [
            (time * MICROSECONDS_TO_SECONDS, frame, update)
            for time, frame, update in entries
        ]

    def reset(self):
        self.next_entry_index = 0
        self.frame_index = 0
        self.time = 0

        # clear simulation
        self.emit(frame=FrameData(), update=None)

    def advance_by_one_step(self):
        self.advance_to_next_entry()

    def advance_by_seconds(self, dt: float):
        next_time = self.time + dt

        while self.entries[self.next_entry_index][0] <= next_time:
            self.advance_to_next_entry()

        # stall one second before looping to beginning
        last_time = self.entries[-1][0]
        self.time = next_time % (last_time + 1)

    def advance_to_next_entry(self):
        time, frame, update = self.entries[self.next_entry_index]
        self.next_entry_index = (self.next_entry_index + 1) % len(self.entries)
        self.time = time
        self.emit(frame=frame, update=update)

    def emit(self, *, frame: Optional[FrameData], update: Optional[DictionaryChange]):
        if self.app_server is None:
            return

        if frame is not None:
            self.app_server.frame_publisher.send_frame(self.frame_index, frame)
            self.frame_index += 1
        if update is not None:
            self.app_server.server.update_state(None, update)
