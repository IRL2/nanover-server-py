"""
Command line interface for connecting to a narupa server and osc server and
generating outgoing osc messages from incoming narupa frame data.
"""
import argparse
import textwrap

from osc_client import OscClient


def frame_to_osc_messages(frame):
    """
    Take a Narupa frame and generate a number of address, message pairs to be
    sent over OSC.
    :param frame: a Narupa frame
    :return: An iterator over OSC address, message pairs
    """
    if frame is None:
        return

    yield "/filter", frame.particle_positions[0]


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
    parser.add_argument('-o', '--osc-port', type=int, default=None)
    parser.add_argument('-i', '--send-interval', type=float, default=.01)
    arguments = parser.parse_args(args)
    return arguments


def initialise(args=None):
    arguments = handle_user_arguments(args)

    runner = OscClient(traj_address=arguments.traj_address,
                       traj_port=arguments.traj_port,
                       osc_address=arguments.osc_address,
                       osc_port=arguments.osc_port,
                       send_interval=arguments.send_interval,
                       message_generator=frame_to_osc_messages)
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
