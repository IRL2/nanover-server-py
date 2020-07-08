# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.
"""
Command line interface for narupa.openmm.
"""
import time
import textwrap
import argparse
from . import OpenMMRunner


def handle_user_arguments() -> argparse.Namespace:
    """
    Parse the arguments from the command line.

    :return: The namespace of arguments read from the command line.
    """
    description = textwrap.dedent("""\
    Run an OpenMM simulation and send it to the network for Narupa.
    """)
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument(
        'simulation_xml_path',
        help='The simulation to run in XML format.',
    )
    parser.add_argument(
        '-v', '--verbose', type=int, default=0, const=100, nargs='?',
        help=('Display the step number, the potential energy in kJ/mol, '
              'and the performance in ns/day.'),
    )
    parser.add_argument(
        '-n', '--name',
        type=str, default=None,
        help='Give a friendly name to the server.',
    )
    parser.add_argument('-p', '--port', type=int, default=None)
    parser.add_argument('-a', '--address', default=None)
    parser.add_argument(
        '-f', '--frame-interval', type=int, default=5, metavar='STEPS',
        help='Sends a frame every STEPS dynamics steps.',
    )
    parser.add_argument(
        '-i', '--force-interval', type=int, default=10, metavar='STEPS',
        help='Update the interactions every STEPS dynamics steps.',
    )

    arguments = parser.parse_args()
    return arguments


def main():
    """
    Entry point for the command line.
    """
    arguments = handle_user_arguments()

    runner = OpenMMRunner.from_xml_input(
        input_xml=arguments.simulation_xml_path,
        name=arguments.name,
        address=arguments.address,
        port=arguments.port,
    )
    print(
        f'Serving "{runner.app.name}" on port {runner.app.port}, '
        f'discoverable on all interfaces on port {runner.app.discovery.port}'
    )

    with runner:
        runner.verbosity_interval = arguments.verbose
        runner.frame_interval = arguments.frame_interval
        runner.force_interval = arguments.force_interval
        runner.run()

        try:
            while True:
                time.sleep(1)
        except KeyboardInterrupt:
            print("Closing due to keyboard interrupt.")


if __name__ == '__main__':
    main()
