import json
import urllib
from contextlib import contextmanager

from websockets.sync.client import connect

from nanover.app import NanoverImdApplication


class DiscoveryClient:
    @classmethod
    @contextmanager
    def advertise_server(cls, endpoint, *, app_server: NanoverImdApplication):
        assert app_server.running_discovery

        ip = get_local_ip()

        data = {
            "name": app_server._service_hub.name,
        }

        services = app_server._service_hub.properties["services"]

        if "wss" in services:
            data["wss"] = f"wss://{ip}:{services["wss"]}"
        if "ws" in services:
            data["ws"] = f"ws://{ip}:{services["ws"]}"

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


def get_local_ip():
    import socket

    def attempt():
        yield from (
            ip
            for ip in socket.gethostbyname_ex(socket.gethostname())[2]
            if ip.startswith("192.")
        )

        with socket.socket(socket.AF_INET, socket.SOCK_DGRAM) as s:
            s.connect(("8.8.8.8", 53))
            name = s.getsockname()[0]
            s.close()
        yield name

    ip = next(attempt())

    return ip
