from nanover.core import AppServerMinimalImd
from nanover.trajectory import FrameData

from .frame_listener import FrameListener
from .utilities import InteractionsUtility


class ImdAgent(FrameListener):
    """
    Base class for agent that attaches to the iMD runner and updates interactions live as simulation frames are
    published.
    """

    def __init__(self, app_server: AppServerMinimalImd):
        super().__init__(app_server)
        self.interactions = InteractionsUtility(app_server)
        self.paused = False

    def on_frame_update(self, full_frame: FrameData, frame_update: FrameData):
        if not self.paused:
            self.update_interactions(full_frame=full_frame, frame_update=frame_update)
        else:
            self.interactions.clear()

    def update_interactions(self, full_frame: FrameData, frame_update: FrameData):
        """
        Update this agent's interactions based on a given frame.
        """

    def close(self):
        """
        Cancel subscription to frame updates and remove all interactions.
        """
        super().close()
        self.interactions.clear()
