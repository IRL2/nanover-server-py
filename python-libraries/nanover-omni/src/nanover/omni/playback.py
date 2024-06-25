from pathlib import Path
from queue import Queue
from typing import List, Optional

from nanover.app import NanoverImdApplication
from nanover.mdanalysis.recordings import iter_trajectory_file, FrameEntry
from nanover.utilities.timing import yield_interval


class PlaybackSimulation:
    def __init__(self, paths: List[str]):
        paths = [Path(path) for path in paths]

        self.name = paths[0].stem

        self.traj_path = next((path for path in paths if path.suffix == ".traj"), None)
        self.state_path = next(
            (path for path in paths if path.suffix == ".state"), None
        )

        self.app_server: Optional[NanoverImdApplication] = None

        self.frames: List[FrameEntry] = []
        self.frame_index = 0
        self.time = 0

    def load(self):
        self.frames = []
        self.frame_index = 0
        self.time = 0

        if self.traj_path:
            self.frames = [
                (elapsed / 1000000, index, frame)
                for elapsed, index, frame in iter_trajectory_file(self.traj_path)
            ]

        if self.state_path:
            pass

    def run(self, app_server: NanoverImdApplication, cancel: Queue):
        self.load()

        last_time = self.frames[-1][0]

        for dt in yield_interval(1 / 30):
            if not cancel.empty():
                break

            prev_time = self.time
            next_time = prev_time + dt

            for time, index, frame in self.frames:
                if prev_time <= time < next_time:
                    app_server.frame_publisher.send_frame(self.frame_index, frame)
                    self.frame_index += 1

            # loop one second after the last frame
            self.time = next_time % (last_time + 1)
