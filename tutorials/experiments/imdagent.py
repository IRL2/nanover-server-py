from concurrent.futures import ThreadPoolExecutor

from nanover.app import OmniRunner
from nanover.imd import ParticleInteraction
from nanover.trajectory import FrameData
from nanover.utilities.cli import CancellationToken


class ImdAgent:
    @classmethod
    def from_runner(cls, runner: OmniRunner):
        return cls(runner)

    def __init__(self, runner: OmniRunner):
        self._runner = runner
        self._threads = ThreadPoolExecutor(max_workers=1)
        self._cancellation = CancellationToken()
        self._task = None
        self._interactions: set[str] = set()

    def update_interaction(self, key: str, interaction: ParticleInteraction):
        self._runner.app_server.imd.insert_interaction(key, interaction)
        self._interactions.add(key)

    def remove_interaction(self, key):
        self._runner.app_server.imd.remove_interaction(key)
        self._interactions.remove(key)

    def update_interactions(self, frame: FrameData):
        pass

    def start(self):
        if self._task is not None:
            return

        publisher = self._runner.app_server.frame_publisher
        stream = publisher.subscribe_latest_frames(
            frame_interval=0, cancellation=self._cancellation
        )

        def run():
            for frame in stream:
                self.update_interactions(frame)

        self._task = self._threads.submit(run)

    def stop(self):
        self._cancellation.cancel()
        self._threads.shutdown(wait=True)
        for key in self._interactions:
            self._runner.app_server.imd.remove_interaction(key)