from narupa.app import NarupaClient
from narupa.core.timing import yield_interval
from pythonosc import udp_client

DEFAULT_OSC_ADDRESS = 'localhost'
DEFAULT_OSC_PORT = 60000


def null_message_generator(frame):
    pass


class OscClient:
    def __init__(self,
                 *,
                 osc_address=None, osc_port=None,
                 traj_address=None, traj_port=None,
                 send_interval=1/30,
                 message_generator=null_message_generator,
                 verbose=False):
        osc_address = osc_address or DEFAULT_OSC_ADDRESS
        osc_port = osc_port or DEFAULT_OSC_PORT

        self.verbose = verbose
        self.message_generator = message_generator
        self.send_interval = send_interval
        self.osc_client = udp_client.SimpleUDPClient(osc_address,
                                                     osc_port,
                                                     allow_broadcast=True)
        self.frame_client = NarupaClient(address=traj_address,
                                         trajectory_port=traj_port)

    def run(self):
        for dt in yield_interval(self.send_interval):
            frame = self.frame_client.latest_frame
            if frame is not None:
                self.process_frame(frame)

    def close(self):
        self.frame_client.close()

    def process_frame(self, frame):
        for address, message in self.message_generator(frame):
            self.osc_client.send_message(address, message)
            if self.verbose:
                print(address, message)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
