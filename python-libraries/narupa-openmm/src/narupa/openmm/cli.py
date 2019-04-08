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
    parser.add_argument('-p', '--port', default=54321)
    parser.add_argument('-a', '--address', default='[::]')
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
            port=arguments.port,
        )
    else:
        simulation_runner = Runner.from_xml_input(arguments.simulation_xml_path)
    simulation_runner.verbose = arguments.verbose
    simulation_runner.run()


if __name__ == '__main__':
    main()
