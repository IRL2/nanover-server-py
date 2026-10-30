import json
import urllib.request
from dataclasses import dataclass
from typing import TypedDict

from websockets.sync.client import ClientConnection, connect

from nanover.core import AppServer
from nanover.utilities.network import get_local_ip


class DiscoveryEndpoint:
    @classmethod
    def from_host(cls, host: str):
        return cls(host)

    def __init__(self, host: str):
        self.host = host

    def get_listing(self):
        endpoint = f"https://{self.host}/list"
        content = urllib.request.urlopen(endpoint).read()
        return json.loads(content)

    def advertise_server(self, app_server: AppServer):
        ip = get_local_ip()

        data = {"name": app_server.service_hub.name}
        services = app_server.service_hub.properties["services"]

        for protocol in ("wss", "ws", "https"):
            if protocol in services:
                port = services[protocol]
                data[protocol] = f"{protocol}://{ip}:{port}"

        return self.advertise_data(data)

    def advertise_data(self, data: dict):
        return DiscoveryAdvert.advertise(data=data, host=self.host)


class DiscoveryProfile(TypedDict):
    id: str


@dataclass
class DiscoveryAdvert:
    host: str
    data: dict
    profile: DiscoveryProfile
    websocket: ClientConnection

    @classmethod
    def advertise(cls, *, data: dict, host: str):
        websocket = None
        try:
            websocket = connect(f"wss://{host}/")
            profile = json.loads(websocket.recv())
            websocket.send(json.dumps(data))
        except Exception as e:
            if websocket is not None:
                websocket.close()
            raise RuntimeError(f"Failed to advertise to {host}") from e

        return cls(host=host, data=data, profile=profile, websocket=websocket)

    def close(self):
        self.websocket.close()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
