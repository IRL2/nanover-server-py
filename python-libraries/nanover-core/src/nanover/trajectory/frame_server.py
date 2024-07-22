from typing import Optional
from nanover.core import (
    NanoverServer,
    get_requested_port_or_default,
    DEFAULT_SERVE_ADDRESS,
)
from .frame_data import FrameData
from .frame_publisher import FramePublisher

DEFAULT_PORT = 54321

PLAY_COMMAND_KEY = "playback/play"
RESET_COMMAND_KEY = "playback/reset"
STEP_COMMAND_KEY = "playback/step"
STEP_BACK_COMMAND_KEY = "playback/step_back"
PAUSE_COMMAND_KEY = "playback/pause"
GET_DYNAMICS_INTERVAL_COMMAND_KEY = "trajectory/get-dynamics-interval"
SET_DYNAMICS_INTERVAL_COMMAND_KEY = "trajectory/set-dynamics-interval"
LIST_COMMAND_KEY = "playback/list"
LOAD_COMMAND_KEY = "playback/load"
NEXT_COMMAND_KEY = "playback/next"


class FrameServer(NanoverServer):
    _trajectory_service: FramePublisher

    def __init__(self, *, address: Optional[str] = None, port: Optional[int] = None):
        if address is None:
            address = DEFAULT_SERVE_ADDRESS
        port = get_requested_port_or_default(port, DEFAULT_PORT)
        super().__init__(address=address, port=port)
        self._frame_count = 0

    def setup_services(self):
        super().setup_services()
        self._trajectory_service = FramePublisher()
        self._trajectory_service.add_to_server_method(
            self._trajectory_service, self.server
        )

    def send_frame(self, frame_index: int, frame_data: FrameData):
        self._trajectory_service.send_frame(frame_index, frame_data.raw)
        self._frame_count += 1

    def close(self):
        super().close()
        self._trajectory_service.close()

    @property
    def frame_count(self):
        """
        Counts how many times send_frame has been called on this publisher.
        """
        return self._frame_count
