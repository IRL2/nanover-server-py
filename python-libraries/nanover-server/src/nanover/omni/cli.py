"""
Command line interface for nanover.omni.
"""

import logging
import time
import textwrap
import argparse
from contextlib import contextmanager
from glob import glob
from typing import Iterable

from nanover.omni import OmniRunner
from nanover.omni.openmm import OpenMMSimulation
from nanover.omni.playback import PlaybackSimulation
from nanover.omni.record import record_from_server
from nanover.utilities.cli import suppress_keyboard_interrupt_as_cancellation


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
        nargs="+",
        default=[],
        metavar="PATH",
        help="Simulation to run via OpenMM (XML format)",
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
        nargs="?",
        default=None,
        const="",
        metavar="PATH",
        help="Record trajectory and state to files.",
    )

    parser.add_argument(
        "-q",
        "--include-velocities",
        action="store_true",
        default=False,
        help="Optionally include the particle velocities in the frame data for OMM and ASE simulations.",
    )
    parser.add_argument(
        "-k",
        "--include-forces",
        action="store_true",
        default=False,
        help="Optionally include the particle forces in the frame data for OMM and ASE simulations.",
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


def get_all_paths(path_sets: Iterable[Iterable[str]]):
    for paths in path_sets:
        for pattern in paths:
            files = list(glob(pattern, recursive=True))
            if files:
                yield from files
            else:
                logging.warning(f'Path "{pattern}" yielded 0 files.')


@contextmanager
def initialise_runner(arguments: argparse.Namespace):
    with OmniRunner.with_basic_server(
        name=arguments.name,
        address=arguments.address,
        port=arguments.port,
    ) as runner:
        for paths in arguments.recording_entries:
            runner.add_simulation(PlaybackSimulation.from_paths(paths))

        for path in get_all_paths(arguments.openmm_xml_entries):
            simulation = OpenMMSimulation.from_xml_path(path)
            simulation.include_velocities = arguments.include_velocities
            simulation.include_forces = arguments.include_forces
            runner.add_simulation(simulation)

        if arguments.record_to_path is not None:
            stem = arguments.record_to_path
            if stem == "":
                timestamp = time.strftime("%Y-%m-%d-%H%M-%S", time.gmtime())
                stem = f"omni-recording-{timestamp}"

            traj_path = f"{stem}.traj"
            state_path = f"{stem}.state"
            print(f"Recording to {traj_path} & {state_path}")

            record_from_server(
                f"localhost:{runner.app_server.port}",
                traj_path,
                state_path,
            )

        yield runner


def main():
    """
    Entry point for the command line.
    """
    # use the nice logger formatting if available
    try:
        from rich.logging import RichHandler

        logging.basicConfig(handlers=[RichHandler(rich_tracebacks=True)])
    except ImportError:
        logging.basicConfig()
    logging.captureWarnings(True)

    arguments = handle_user_arguments()

    with suppress_keyboard_interrupt_as_cancellation() as cancellation:
        with initialise_runner(arguments) as runner:
            if len(runner.simulations) > 0:
                runner.load(0)

            runner.print_basic_info()
            cancellation.wait_cancellation(interval=0.5)
            print("Closing due to KeyboardInterrupt.")


if __name__ == "__main__":
    main()
