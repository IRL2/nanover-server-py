# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.
"""
Command line interface for narupa.openmm.
"""
import textwrap
import argparse
from . import Runner


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
        '-v', '--verbose',
        action='store_true',
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

    simulation_runner = Runner.from_xml_input(
        input_xml=arguments.simulation_xml_path,
        name=arguments.name,
        address=arguments.address,
        port=arguments.port,
    )
    print(f'Serving frames on port {simulation_runner.app.port}.')

    with simulation_runner:
        simulation_runner.verbose = arguments.verbose
        simulation_runner.frame_interval = arguments.frame_interval
        simulation_runner.force_interval = arguments.force_interval
        try:
            simulation_runner.run()
        except KeyboardInterrupt:
            print("Closing due to keyboard interrupt.")


if __name__ == '__main__':
    main()
