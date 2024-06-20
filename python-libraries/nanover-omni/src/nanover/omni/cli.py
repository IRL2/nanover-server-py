"""
Command line interface for nanover.omni.
"""

import time
import textwrap
import argparse

from nanover.app import NanoverImdApplication


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
        "openmm_xml_paths"
        "--omm",
        "--openmm",
        action="append",
        nargs=1,
        help="Simulation to run via OpenMM (XML format)",
    )

    parser.add_argument(
        "ase_xml_paths"
        "--ase",
        "--ase-omm",
        "--ase-openmm",
        action="append",
        nargs=1,
        help="Simulation to run via ASE OpenMM (XML format)",
    )

    parser.add_argument(
        ""
        "--playback",
        "--recording",
        action="append",
        nargs="+",
        help="Recorded session to playback (one or both of .traj and .state)",
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

    arguments = parser.parse_args()
    return arguments


def main():
    """
    Entry point for the command line.
    """
    arguments = handle_user_arguments()

    print(arguments)

    app_server = NanoverImdApplication.basic_server(
        name=arguments.name,
        address=arguments.address,
        port=arguments.port,
    )

    return


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
