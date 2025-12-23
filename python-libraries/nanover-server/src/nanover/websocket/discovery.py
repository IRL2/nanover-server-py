import json
import urllib
from contextlib import contextmanager

from websockets.sync.client import connect

from nanover.core import AppServer
from nanover.utilities.network import get_local_ip


class DiscoveryClient:
    @classmethod
    @contextmanager
    def advertise_server(cls, endpoint, *, app_server: AppServer):
        ip = get_local_ip()

        data = {
            "name": app_server.service_hub.name,
        }

        services = app_server.service_hub.properties["services"]

        if "wss" in services:
            data["wss"] = f"wss://{ip}:{services['wss']}"
        if "ws" in services:
            data["ws"] = f"ws://{ip}:{services['ws']}"

        discovery = DiscoveryClient(endpoint)
        with discovery.advertise(data) as init:
            yield discovery, init

    def __init__(self, host: str):
        self.host = host

    def get_listing(self):
        endpoint = f"https://{self.host}/list"
        content = urllib.request.urlopen(endpoint).read()
        return json.loads(content)

    @contextmanager
    def advertise(self, data):
        with connect(f"wss://{self.host}/") as websocket:
            init = json.loads(websocket.recv())
            websocket.send(json.dumps(data))
            yield init
