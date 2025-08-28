"""
Command line interface for nanover.omni.
"""

import logging
import ssl
import time
import textwrap
import argparse
from contextlib import contextmanager, nullcontext
from glob import glob
from typing import Iterable

from nanover.essd import ServiceHub
from nanover.omni import OmniRunner
from nanover.omni.openmm import OpenMMSimulation
from nanover.omni.playback import PlaybackSimulation
from nanover.omni.record import record_from_server
from nanover.utilities.cli import suppress_keyboard_interrupt_as_cancellation
from nanover.websocket.discovery import get_local_ip, DiscoveryClient


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
        ssl=initialise_ssl(arguments),
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


def initialise_ssl(arguments: argparse.Namespace):
    if arguments.ssl is None:
        return None

    certfile, keyfile, password = arguments.ssl
    ssl_context = ssl.SSLContext(ssl.PROTOCOL_TLS_SERVER)
    ssl_context.load_verify_locations(certfile)
    ssl_context.load_cert_chain(certfile, keyfile=keyfile, password=password)

    return ssl_context


@contextmanager
def do_cloud_discovery(endpoint, *, hub: ServiceHub):
    ip = get_local_ip()

    data = {
        "name": hub.name,
        "web": f"https://{ip}:5500",
        "https": f"https://{ip}:5500",
    }

    services = hub.properties["services"]

    if "wss" in services:
        data["wss"] = f"wss://{ip}:{services["wss"]}"
    if "ws" in services:
        data["ws"] = f"ws://{ip}:{services["ws"]}"

    print(data)
    discovery = DiscoveryClient(endpoint)
    print(discovery.get_listing())
    with discovery.advertise(data) as init:
        yield init


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

            if arguments.cloud_discovery_host is not None:
                cloud_discovery = do_cloud_discovery(
                    arguments.cloud_discovery_host, hub=runner.app_server._service_hub
                )
            else:
                cloud_discovery = nullcontext()  # type: ignore

            with cloud_discovery:
                cancellation.wait_cancellation()
            print("Closing due to KeyboardInterrupt.")


if __name__ == "__main__":
    main()
