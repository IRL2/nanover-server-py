from concurrent.futures import Future, ThreadPoolExecutor
from logging import getLogger

from nanover.app import OmniRunner
from nanover.core import AppServerMinimalImd
from nanover.trajectory import FrameData
from nanover.utilities.cli import CancellationToken


class FrameListener:
    """
    Base class for agent that attaches to the iMD runner and calls an update method when frames are
    published.
    """

    @classmethod
    def from_runner(cls, runner: OmniRunner):
        return cls(runner.app_server)

    def __init__(self, app_server: AppServerMinimalImd):
        self._app_server = app_server
        self._threads = ThreadPoolExecutor(max_workers=1)
        self._cancellation = CancellationToken()
        self._task: Future | None = None

    def on_frame_update(self, full_frame: FrameData, frame_update: FrameData):
        pass

    def start(self):
        """
        Subscribe to frame updates and update interactions until closed or paused.
        """
        if self._task is not None:
            return

        publisher = self._app_server.frame_publisher
        stream = publisher.subscribe_latest_frames(
            frame_interval=0, cancellation=self._cancellation
        )

        def run():
            full_frame = FrameData()
            for frame_update in stream:
                full_frame.update(frame_update)
                try:
                    self.on_frame_update(
                        full_frame=full_frame,
                        frame_update=frame_update,
                    )
                except Exception:
                    getLogger().exception("Exception in `on_frame_update`")

        self._task = self._threads.submit(run)

    def close(self):
        """
        Cancel subscription to frame.
        """
        self._cancellation.cancel()
        self._threads.shutdown(wait=True)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
