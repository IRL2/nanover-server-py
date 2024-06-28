"""
Command line interface for nanover.omni.
"""

import time
import textwrap
import argparse

from nanover.app import NanoverImdApplication
from nanover.omni import OmniRunner
from nanover.omni.openmm import OpenMMSimulation
from nanover.omni.ase_omm import ASEOpenMMSimulation
from nanover.omni.playback import PlaybackSimulation


def handle_user_arguments() -> argparse.Namespace:
    """
    Parse the arguments from the command line.

    :return: The namespace of arguments read from the command line.
    """
    description = textwrap.dedent(
        """\
    Multi-simulation server for NanoVer
    """
    )
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument(
        "--omm",
        "--openmm",
        dest="openmm_xml_entries",
        action="append",
        nargs=1,
        default=[],
        help="Simulation to run via OpenMM (XML format)",
    )

    parser.add_argument(
        "--ase",
        "--ase-omm",
        "--ase-openmm",
        dest="ase_xml_entries",
        action="append",
        nargs=1,
        default=[],
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

    for paths in arguments.recording_entries:
        runner.add_simulation(PlaybackSimulation.from_paths(paths))

    for (path,) in arguments.openmm_xml_entries:
        runner.add_simulation(OpenMMSimulation(path))

    for (path,) in arguments.ase_xml_entries:
        runner.add_simulation(ASEOpenMMSimulation(path))

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
        runner.next()
        try:
            while True:
                time.sleep(1)
        except KeyboardInterrupt:
            print("Closing due to keyboard interrupt.")


if __name__ == "__main__":
    main()
