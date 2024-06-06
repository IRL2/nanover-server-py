"""
Command line interface for nanover.openmm.
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
    description = textwrap.dedent(
        """\
    Run an OpenMM simulation and send it to the network for NanoVer.
    """
    )
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument(
        "simulation_xml_paths",
        nargs="+",
        help="The simulations to run in XML format.",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        type=int,
        default=0,
        const=100,
        nargs="?",
        help=(
            "Display the step number, the potential energy in kJ/mol, "
            "and the performance in ns/day."
        ),
    )
    parser.add_argument(
        "-n",
        "--name",
        type=str,
        default=None,
        help="Give a friendly name to the server.",
    )
    parser.add_argument("-p", "--port", type=int, default=None)
    parser.add_argument("-a", "--address", default=None)
    parser.add_argument(
        "-f",
        "--frame-interval",
        type=int,
        default=5,
        metavar="STEPS",
        help="Sends a frame every STEPS dynamics steps.",
    )
    parser.add_argument(
        "-i",
        "--force-interval",
        type=int,
        default=10,
        metavar="STEPS",
        help="Update the interactions every STEPS dynamics steps.",
    )
    parser.add_argument(
        "--platform", default=None, help="Select the platform on which to run Openmm."
    )

    arguments = parser.parse_args()
    return arguments


def main():
    """
    Entry point for the command line.
    """
    arguments = handle_user_arguments()

    runner = OpenMMRunner.from_xml_inputs(
        input_xmls=arguments.simulation_xml_paths,
        name=arguments.name,
        address=arguments.address,
        port=arguments.port,
        platform=arguments.platform,
    )
    print(
        f'Serving "{runner.app_server.name}" on port {runner.app_server.port}, '
        f"discoverable on all interfaces on port {runner.app_server.discovery.port}"
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


if __name__ == "__main__":
    main()
