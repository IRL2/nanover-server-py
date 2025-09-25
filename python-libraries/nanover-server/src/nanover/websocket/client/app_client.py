import time
from collections import deque
from typing import Any

from nanover.app.types import AppServer
from nanover.essd import DiscoveryClient, ServiceHub
from nanover.utilities.change_buffers import DictionaryChange
from nanover.websocket.client.playback_client import PlaybackClient
from nanover.utilities.network import get_local_ip
from nanover.websocket.client.interaction_client import InteractionClient
from nanover.websocket.client.selection_client import SelectionClient
from nanover.trajectory import FrameData2

DEFAULT_DISCOVERY_SEARCH_TIME = 10.0


class NanoverImdClient(InteractionClient, SelectionClient, PlaybackClient):
    @classmethod
    def from_runner(cls, runner: Any):
        return cls.from_app_server(runner.app_server)

    @classmethod
    def from_app_server(cls, app_server: AppServer):
        try:
            return cls.from_url(
                f"wss://localhost:{app_server.service_hub.services["wss"]}"
            )
        except KeyError:
            return cls.from_url(
                f"ws://localhost:{app_server.service_hub.services["ws"]}"
            )

    @classmethod
    def from_discovery(
        cls,
        *,
        server_name: str | None = None,
        discovery_port: int | None = None,
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
        self._frames: deque[FrameData2] = deque(maxlen=50)
        super().__init__(*args, **kwargs)

    @property
    def frames(self) -> list[FrameData2]:
        return list(self._frames)

    def recv_frame(self, message: dict):
        super().recv_frame(message)
        self._frames.append(FrameData2(message).copy())

    @property
    def current_frame(self) -> FrameData2:
        return self._current_frame

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


def get_websocket_address_from_hub(hub: ServiceHub, *, host: str | None = None):
    parts = hub.get_service_address("ws")
    assert parts is not None
    host_, port = parts
    address = f"ws://{host or host_}:{port}"
    return address


def get_websocket_address_from_app_server(app_server: AppServer):
    return get_websocket_address_from_hub(app_server.service_hub, host="localhost")


def _search_for_first_server_with_name(
    server_name: str,
    search_time: float = DEFAULT_DISCOVERY_SEARCH_TIME,
    discovery_address: str | None = None,
    discovery_port: int | None = None,
):
    with DiscoveryClient(discovery_address, discovery_port) as discovery_client:
        for hub in discovery_client.search_for_services(search_time):
            if hub.name == server_name:
                return hub
    return None


def _search_for_first_local_server(
    search_time: float = DEFAULT_DISCOVERY_SEARCH_TIME,
    discovery_address: str | None = None,
    discovery_port: int | None = None,
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
