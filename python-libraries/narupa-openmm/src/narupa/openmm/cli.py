"""
Command line interface for narupa.openmm.
"""
import argparse
from . import Runner


def handle_user_arguments() -> argparse.Namespace:
    """
    Parse the arguments from the command line.

    :return: The namespace of arguments read from the command line.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('system_xml_path')
    parser.add_argument('-v', '--verbose', action='store_true')

    arguments = parser.parse_args()
    return arguments


def main():
    """
    Entry point for the command line.
    """
    arguments = handle_user_arguments()
    simulation_runner = Runner.from_xml_input(
        input_xml=arguments.system_xml_path,
    )
    simulation_runner.verbose = arguments.verbose
    simulation_runner.run()


if __name__ == '__main__':
    main()
