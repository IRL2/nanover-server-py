from concurrent.futures import ThreadPoolExecutor, Future

from nanover.app import OmniRunner
from nanover.core import AppServerMinimalImd
from nanover.imd import ParticleInteraction, INTERACTION_PREFIX
from nanover.trajectory import FrameData
from nanover.utilities.cli import CancellationToken
from nanover.websocket.client.app_client import NanoverImdClient


class ImdAgent:
    """
    Base class for agent that attaches to the iMD runner and updates interactions live as simulation frames are
    published.
    """

    @classmethod
    def from_client(cls, client: NanoverImdClient):
        return cls.from_app_server(client)

    @classmethod
    def from_runner(cls, runner: OmniRunner):
        return cls(runner.app_server)

    @classmethod
    def from_app_server(cls, app_server: AppServerMinimalImd):
        return cls(app_server)

    def __init__(self, app_server: AppServerMinimalImd):
        self._app_server = app_server
        self._threads = ThreadPoolExecutor(max_workers=1)
        self._cancellation = CancellationToken()
        self._task: Future | None = None
        self._interactions: set[str] = set()
        self.paused = False

    def update_interaction(self, name: str, interaction: ParticleInteraction):
        """
        Add or update a named interaction.
        """
        key = f"{INTERACTION_PREFIX}.{name}"
        self._app_server.imd.insert_interaction(key, interaction)
        self._interactions.add(key)

    def remove_interaction(self, name):
        """
        Remove a named interaction.
        """
        key = f"{INTERACTION_PREFIX}.{name}"
        self._app_server.imd.remove_interaction(key)
        self._interactions.remove(key)

    def clear_interactions(self):
        """
        Clear all interactions created by this agent.
        """
        for key in self._interactions:
            self._app_server.imd.remove_interaction(key)
        self._interactions.clear()

    def update_interactions(self, full_frame: FrameData, frame_update: FrameData):
        """
        Update this agent's interactions based on a given frame.
        """
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
                if not self.paused:
                    self.update_interactions(full_frame, frame_update)
                else:
                    self.clear_interactions()

        self._task = self._threads.submit(run)

    def close(self):
        """
        Cancel subscription to frame updates and remove all interactions.
        """
        self._cancellation.cancel()
        self._threads.shutdown(wait=True)
        self.clear_interactions()
