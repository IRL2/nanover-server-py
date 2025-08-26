import json
import urllib
from contextlib import contextmanager

from websockets.sync.client import connect


class DiscoveryClient:
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
