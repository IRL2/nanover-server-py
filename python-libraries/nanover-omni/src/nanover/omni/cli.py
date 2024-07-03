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
from nanover.omni.record import record_from_server


def handle_user_arguments(args=None) -> argparse.Namespace:
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
        dest="openmm_xml_entries",
        action="append",
        default=[],
        metavar="PATH",
        help="Simulation to run via OpenMM (XML format)",
    )

    parser.add_argument(
        "--ase-omm",
        dest="ase_xml_entries",
        action="append",
        default=[],
        metavar="PATH",
        help="Simulation to run via ASE OpenMM (XML format)",
    )

    parser.add_argument(
        "--playback",
        dest="recording_entries",
        action="append",
        nargs="+",
        default=[],
        metavar="PATH",
        help="Recorded session to playback (one or both of .traj and .state)",
    )

    parser.add_argument(
        "--record",
        dest="record_to_path",
        metavar="PATH",
        help="Filename to record trajectory and state updates to.",
    )

    parser.add_argument(
        "-n",
        "--name",
        help="Give a friendly name to the server.",
    )
    parser.add_argument("-p", "--port", type=int, default=None)
    parser.add_argument("-a", "--address", default=None)

    arguments = parser.parse_args(args)
    return arguments


def initialise(args=None):
    arguments = handle_user_arguments(args)

    runner = OmniRunner.with_basic_server(
        name=arguments.name,
        address=arguments.address,
        port=arguments.port,
    )

    for paths in arguments.recording_entries:
        runner.add_simulation(PlaybackSimulation.from_paths(paths))

    for path in arguments.openmm_xml_entries:
        runner.add_simulation(OpenMMSimulation(path))

    for path in arguments.ase_xml_entries:
        runner.add_simulation(ASEOpenMMSimulation(path))

    if arguments.record_to_path is not None:
        traj_path = f"{arguments.record_to_path}.traj"
        state_path = f"{arguments.record_to_path}.state"
        print(f"Recording to {traj_path} & {state_path}")

        record_from_server(
            f"localhost:{runner.app_server.port}",
            traj_path,
            state_path,
        )

    return runner


def main():
    """
    Entry point for the command line.
    """
    with initialise() as runner:
        print(
            f'Serving "{runner.app_server.name}" on port {runner.app_server.port}, '
            f"discoverable on all interfaces on port {runner.app_server.discovery.port}"
        )

        list = "\n".join(
            f'{index}: "{simulation.name}"'
            for index, simulation in enumerate(runner.simulations)
        )
        print(f"Available simulations:\n{list}")

        if len(runner.simulations) > 0:
            runner.next()
        try:
            while True:
                time.sleep(1)
        except KeyboardInterrupt:
            print("Closing due to keyboard interrupt.")


if __name__ == "__main__":
    main()
