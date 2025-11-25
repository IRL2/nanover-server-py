"""
Command line interface for nanover.omni.
"""

import logging
import ssl
import time
import textwrap
import argparse
from contextlib import contextmanager
from glob import glob
from typing import Iterable

from nanover.omni import OmniRunner
from nanover.mdanalysis import UniverseSimulation
from nanover.openmm import OpenMMSimulation
from nanover.recording import PlaybackSimulation
from nanover.utilities.cli import suppress_keyboard_interrupt_as_cancellation
from nanover.websocket.client.app_client import NanoverImdClient
from nanover.websocket.discovery import DiscoveryClient
from nanover.websocket.record import record_from_runner


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
        help="Simulations to run via OpenMM (XML format)",
    )

    parser.add_argument(
        "--mda",
        dest="mdanalysis_entries",
        action="append",
        nargs="+",
        default=[],
        metavar="PATH",
        help="Structures to load via MDanalysis",
    )

    parser.add_argument(
        "--playback",
        dest="recording_entries",
        action="append",
        nargs="+",
        default=[],
        metavar="PATH",
        help="Recorded sessions to playback",
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
        "-s",
        "--ssl",
        nargs=3,
        metavar=("CERTFILE", "KEYFILE", "PASSWORD"),
        help="Offer SSL in addition to insecure connections.",
    )

    parser.add_argument(
        "--cloud-discovery",
        dest="cloud_discovery_host",
        metavar="ENDPOINT",
        help="Advertise server on cloud discovery.",
    )

    parser.add_argument(
        "--remote",
        nargs="?",
        metavar="PROTOCOL://HOST:PORT",
        help="Publish to remote server at specific host address and port.",
        const="ws://localhost:38801",
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
    @contextmanager
    def make_runner():
        if arguments.remote is None:
            with OmniRunner.with_basic_server(
                name=arguments.name,
                address=arguments.address,
                port=arguments.port,
                ssl=initialise_ssl(arguments),
            ) as runner:
                yield runner
        else:
            with NanoverImdClient.from_url(arguments.remote) as client:
                print(f"Serving as client of remote server at {arguments.remote}")
                with OmniRunner.from_client(client) as runner:
                    yield runner

    with make_runner() as runner:
        for path in get_all_paths(arguments.recording_entries):
            runner.add_simulation(PlaybackSimulation.from_path(path=path))

        for path in get_all_paths(arguments.openmm_xml_entries):
            simulation = OpenMMSimulation.from_xml_path(path)
            simulation.include_velocities = arguments.include_velocities
            simulation.include_forces = arguments.include_forces
            runner.add_simulation(simulation)

        for path in get_all_paths(arguments.mdanalysis_entries):
            simulation = UniverseSimulation.from_path(path=path)
            runner.add_simulation(simulation)

        if arguments.record_to_path is not None:
            stem = arguments.record_to_path
            if stem == "":
                timestamp = time.strftime("%Y-%m-%d-%H%M-%S", time.gmtime())
                stem = f"omni-recording-{timestamp}"

            out_path = f"{stem}.nanover.zip"
            print(f"Recording to {out_path}")

            record_from_runner(runner, out_path)

        if arguments.cloud_discovery_host:
            with DiscoveryClient.advertise_server(
                arguments.cloud_discovery_host, app_server=runner.app_server
            ):
                yield runner
        else:
            yield runner


def initialise_ssl(arguments: argparse.Namespace):
    if arguments.ssl is None:
        return None

    certfile, keyfile, password = arguments.ssl
    ssl_context = ssl.SSLContext(ssl.PROTOCOL_TLS_SERVER)
    ssl_context.load_verify_locations(certfile)
    ssl_context.load_cert_chain(certfile, keyfile=keyfile, password=password)

    return ssl_context


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
            cancellation.wait_cancellation()
            print("Closing due to KeyboardInterrupt.")


def deprecated():
    import warnings

    warnings.warn(
        "use `nanover-server` instead",
        DeprecationWarning,
        stacklevel=2,
    )
    main()


if __name__ == "__main__":
    main()
