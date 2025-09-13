import time
from collections import deque
from typing import Optional, Any

from nanover.app import NanoverImdApplication
from nanover.app.client import _search_for_first_server_with_name
from nanover.essd import DiscoveryClient
from nanover.trajectory import FrameData
from nanover.utilities.change_buffers import DictionaryChange
from nanover.websocket.client.playback_client import PlaybackClient
from nanover.websocket.convert import (
    convert_dict_frame_to_grpc_frame,
    unpack_dict_frame,
)
from nanover.websocket.discovery import get_local_ip
from nanover.websocket.client.interaction_client import InteractionClient
from nanover.websocket.client.selection_client import SelectionClient


class NanoverImdClient(InteractionClient, SelectionClient, PlaybackClient):
    @classmethod
    def from_runner(cls, runner: Any):
        return cls.from_app_server(runner.app_server)

    @classmethod
    def from_app_server(cls, app_server: NanoverImdApplication):
        server = app_server._server_ws
        assert server is not None
        if server.wss_port is not None:
            return cls.from_url(f"wss://localhost:{server.wss_port}")
        else:
            return cls.from_url(f"ws://localhost:{server.ws_port}")

    @classmethod
    def from_discovery(
        cls,
        *,
        server_name: Optional[str] = None,
        discovery_port: Optional[int] = None,
    ):
        if server_name is not None:
            first_service = _search_for_first_server_with_name(
                server_name=server_name, discovery_port=discovery_port
            )
            if first_service is None:
                raise ConnectionError(
                    f'Couldn\'t discover server with name "{server_name}"'
                )
        else:
            first_service = _search_for_first_local_server(
                discovery_port=discovery_port
            )
            if first_service is None:
                raise ConnectionError("Couldn't discover server")

        def get_url(protocol: str):
            address = first_service.get_service_address(service_name=protocol)
            if address is None:
                return None
            hostname, port = address
            return f"{protocol}://{hostname}:{port}"

        url = get_url("wss") or get_url("ws")
        if url is None:
            raise ConnectionError(
                f'Service "{first_service.name}" doesn\'t support wss or ws'
            )

        return cls.from_url(url)

    def __init__(self, *args, **kwargs):
        self._frames_grpc: deque[FrameData] = deque(maxlen=50)
        super().__init__(*args, **kwargs)

    @property
    def frames(self):
        return list(self._frames_grpc)

    def recv_frame(self, message: dict):
        # TODO: don't do this
        frame_grpc = convert_dict_frame_to_grpc_frame(unpack_dict_frame(message))
        self._frames_grpc.append(frame_grpc)
        super().recv_frame(message)

    @property
    def current_frame(self):
        return self._current_frame

    @property
    def current_frame_grpc(self):
        return convert_dict_frame_to_grpc_frame(self.current_frame)

    @property
    def latest_multiplayer_values(self):
        return self._state_dictionary.copy_content()

    def wait_until_first_frame(self, check_interval=0.01, timeout=1):
        """
        Wait until the first frame is received from the server.

        :param check_interval: Interval at which to check if a frame has been
            received.
        :param timeout: Timeout after which to stop waiting for a frame.
        :return: The first :class:`FrameData` received.
        :raises Exception: if no frame is received.
        """
        endtime = 0 if timeout is None else time.monotonic() + timeout

        while not self.current_frame:
            if 0 < endtime < time.monotonic():
                raise Exception("Timed out waiting for first frame.")
            time.sleep(check_interval)

        return self.current_frame

    def update_available_commands(self):
        return self.run_command_blocking("commands/list")["list"]

    def copy_state(self):
        return self._state_dictionary.copy_content()

    def set_shared_value(self, key: str, value):
        self.update_state(DictionaryChange(updates={key: value}))

    def remove_shared_value(self, key: str):
        self.update_state(DictionaryChange(removals={key}))

    def get_shared_value(self, key: str, default=None):
        return self._state_dictionary.copy_content().get(key, default)


def _search_for_first_local_server(
    search_time: float = 2.0,
    discovery_address: Optional[str] = None,
    discovery_port: Optional[int] = None,
):
    try:
        local_ip = get_local_ip()
    except OSError:
        local_ip = None

    with DiscoveryClient(discovery_address, discovery_port) as discovery_client:
        for hub in discovery_client.search_for_services(search_time):
            if hub.address == "localhost" or hub.address == local_ip:
                return hub
    return None
