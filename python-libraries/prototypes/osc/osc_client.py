"""
Command line interface for connecting to a narupa server and osc server and
generating outgoing osc messages from incoming narupa frame data.
"""
import argparse
import textwrap
import time

from narupa.app import NarupaClient
from pythonosc import udp_client

DEFAULT_OSC_ADDRESS = 'localhost'
DEFAULT_OSC_PORT = 60000


def frame_to_osc_messages(frame):
    if frame is None:
        return

    yield "/filter", frame.particle_positions[0]


class OscRunner:
    def __init__(self, args):
        osc_address = args.osc_address or DEFAULT_OSC_ADDRESS
        osc_port = args.osc_port or DEFAULT_OSC_PORT

        self.osc_client = udp_client.SimpleUDPClient(osc_address, osc_port)
        self.frame_client = NarupaClient(address=args.traj_address,
                                         trajectory_port=args.traj_port)

    def run(self):
        while True:
            frame = self.frame_client.latest_frame
            for address, message in frame_to_osc_messages(frame):
                self.osc_client.send_message(address, message)
            time.sleep(1 / 30)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.frame_client.close()


def handle_user_arguments(args=None) -> argparse.Namespace:
    """
    Parse the arguments from the command line.

    :return: The namespace of arguments read from the command line.
    """
    description = textwrap.dedent("""\
    Connect to a narupa trajectory service, and osc server and generate outgoing
    osc messages from incoming frame data.
    """)
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('--traj-address', default=None)
    parser.add_argument('--osc-address', default='localhost')
    parser.add_argument('-t', '--traj-port', type=int, default=None)
    parser.add_argument('-o', '--osc-port', type=int, default=60000)
    arguments = parser.parse_args(args)
    return arguments


def initialise(args=None):
    arguments = handle_user_arguments(args)

    runner = OscRunner(arguments)
    return runner


def main():
    """
    Entry point for the command line.
    """
    with initialise() as runner:
        print('Running...')
        try:
            runner.run()
        except KeyboardInterrupt:
            print("Closing due to keyboard interrupt.")


if __name__ == '__main__':
    main()
