# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.
from narupa.app import NarupaImdClient
from narupa.utilities.timing import yield_interval
from pythonosc import udp_client

# doesn't support both IPv4 and IPv6 at once, so we probably want IPv4
# See https://github.com/attwad/python-osc/issues/109
DEFAULT_OSC_ADDRESS = ('127.0.0.1', 60000)


def null_message_generator(frame):
    pass


class OscClient:
    def __init__(
            self,
            narupa_client: NarupaImdClient,
            *,
            osc_address=DEFAULT_OSC_ADDRESS,
            osc_send_interval=1/30,
            message_generator=null_message_generator,
            verbose=False,
    ):
        self.verbose = verbose
        self.message_generator = message_generator
        self.send_interval = osc_send_interval

        host, port = osc_address
        self.osc_client = udp_client.SimpleUDPClient(host, port, allow_broadcast=True)
        self.narupa_client = narupa_client
        self.narupa_client.subscribe_to_all_frames()

    def run(self):
        for dt in yield_interval(self.send_interval):
            frame = self.narupa_client.latest_frame
            if frame is not None:
                self.process_frame(frame)

    def close(self):
        self.narupa_client.close()

    def process_frame(self, frame):
        for address, message in self.message_generator(frame):
            self.osc_client.send_message(address, message)
            if self.verbose:
                print(address, message)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
