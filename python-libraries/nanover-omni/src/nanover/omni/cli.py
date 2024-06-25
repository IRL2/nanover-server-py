"""
Command line interface for nanover.omni.
"""

import time
import textwrap
import argparse

from nanover.app import NanoverImdApplication
from nanover.omni import OmniRunner
from nanover.omni.playback import PlaybackSimulation


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
        "--omm",
        "--openmm",
        dest="openmm_xml_entries",
        action="append",
        nargs=1,
        help="Simulation to run via OpenMM (XML format)",
    )

    parser.add_argument(
        "--ase",
        "--ase-omm",
        "--ase-openmm",
        dest="ase_xml_entries",
        action="append",
        nargs=1,
        help="Simulation to run via ASE OpenMM (XML format)",
    )

    parser.add_argument(
        "--playback",
        "--recording",
        dest="recording_entries",
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

    app_server = NanoverImdApplication.basic_server(
        name=arguments.name,
        address=arguments.address,
        port=arguments.port,
    )

    runner = OmniRunner(app_server)

    for entry in arguments.recording_entries:
        runner.add_simulation(PlaybackSimulation(entry))

    runner.next()

    print(
        f'Serving "{runner.app_server.name}" on port {runner.app_server.port}, '
        f"discoverable on all interfaces on port {runner.app_server.discovery.port}"
    )

    list = ", ".join(
        f'{index}: "{simulation.name}"'
        for index, simulation in enumerate(runner.simulations)
    )
    print(f"Available simulations: {list}")

    with runner:
        # runner.verbosity_interval = arguments.verbose
        # runner.frame_interval = arguments.frame_interval
        # runner.force_interval = arguments.force_interval
        # runner.run()

        try:
            while True:
                time.sleep(5)
                runner.next()
        except KeyboardInterrupt:
            print("Closing due to keyboard interrupt.")


if __name__ == "__main__":
    main()
