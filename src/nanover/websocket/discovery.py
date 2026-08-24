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

        for protocol in ("wss", "ws", "https"):
            if protocol in services:
                port = services[protocol]
                data[protocol] = f"{protocol}://{ip}:{port}"

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
