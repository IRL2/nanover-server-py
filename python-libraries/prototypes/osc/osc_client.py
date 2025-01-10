from nanover.app import NanoverImdClient
from nanover.utilities.timing import yield_interval
from pythonosc import udp_client

# doesn't support both IPv4 and IPv6 at once, so we probably want IPv4
# See https://github.com/attwad/python-osc/issues/109
DEFAULT_OSC_ADDRESS = ("127.0.0.1", 60000)


def null_message_generator(frame):
    pass


class OscClient:
    def __init__(
        self,
        nanover_client: NanoverImdClient,
        *,
        osc_address=DEFAULT_OSC_ADDRESS,
        osc_send_interval=1 / 30,
        message_generator=null_message_generator,
        verbose=False,
    ):
        self.verbose = verbose
        self.message_generator = message_generator
        self.send_interval = osc_send_interval

        host, port = osc_address
        self.osc_client = udp_client.SimpleUDPClient(host, port, allow_broadcast=True)
        self.nanover_client = nanover_client
        self.nanover_client.subscribe_to_all_frames()

    def run(self):
        for dt in yield_interval(self.send_interval):
            if not self.nanover_client.are_frames_subscribed:
                break

            frame = self.nanover_client.latest_frame
            if frame is not None:
                self.process_frame(frame)

    def close(self):
        self.nanover_client.close()

    def process_frame(self, frame):
        for address, message in self.message_generator(frame):
            self.osc_client.send_message(address, message)
            if self.verbose:
                print(address, message)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
