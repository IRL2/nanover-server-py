# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.
"""
Command line interface for narupa.openmm.
"""
import textwrap
import argparse
from . import Runner, Server


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
    parser.add_argument('-t', '--trajectory_port', default=None)
    parser.add_argument('-a', '--address', default=None)
    parser.add_argument(
        '--no-serve',
        dest='do_serve', action='store_false', default=True,
        help='Do not send the trajectory over the network for Narupa.',
    )

    arguments = parser.parse_args()
    return arguments


def main():
    """
    Entry point for the command line.
    """
    arguments = handle_user_arguments()

    if arguments.do_serve:
        simulation_runner = Server.from_xml_input(
            input_xml=arguments.simulation_xml_path,
            address=arguments.address,
            port=arguments.trajectory_port,
        )
        print(f'Serving frames on port {simulation_runner.trajectory_port}.')
    else:
        simulation_runner = Runner.from_xml_input(arguments.simulation_xml_path)
        print(f'Running without serving frames.')

    with simulation_runner:
        simulation_runner.verbose = arguments.verbose
        try:
            simulation_runner.run()
        except KeyboardInterrupt:
            print("Closing due to keyboard interrupt.")


if __name__ == '__main__':
    main()
