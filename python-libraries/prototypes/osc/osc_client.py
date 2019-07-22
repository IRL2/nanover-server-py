import time

from narupa.app import NarupaClient
from pythonosc import udp_client

DEFAULT_OSC_ADDRESS = 'localhost'
DEFAULT_OSC_PORT = 60000


def null_message_generator(frame):
    pass


def yield_interval(interval):
    """
    Yield at a set interval, accounting for the time spent outside of this
    function.
    :param interval: Number of seconds to put between yields
    """
    last_yield = 0
    while True:
        time_since_yield = time.monotonic() - last_yield
        wait_duration = max(0, interval - time_since_yield)
        time.sleep(wait_duration)
        yield
        last_yield = time.monotonic()


class OscClient:
    def __init__(self,
                 *,
                 osc_address=None, osc_port=None,
                 traj_address=None, traj_port=None,
                 send_interval=1/30,
                 message_generator=None):
        osc_address = osc_address or DEFAULT_OSC_ADDRESS
        osc_port = osc_port or DEFAULT_OSC_PORT

        self.message_generator = message_generator or null_message_generator
        self.send_interval = send_interval
        self.osc_client = udp_client.SimpleUDPClient(osc_address,
                                                     osc_port,
                                                     allow_broadcast=True)
        self.frame_client = NarupaClient(address=traj_address,
                                         trajectory_port=traj_port)

    def run(self):
        for _ in yield_interval(self.send_interval):
            frame = self.frame_client.latest_frame
            if frame is not None:
                self.process_frame(frame)

    def close(self):
        self.frame_client.close()

    def process_frame(self, frame):
        for address, message in self.message_generator(frame):
            self.osc_client.send_message(address, message)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
