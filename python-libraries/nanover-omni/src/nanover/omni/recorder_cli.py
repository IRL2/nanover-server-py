"""
Command line interface for recording from a server.
"""
import argparse
import textwrap
import time

from nanover.omni.record import record_from_server

from nanover.essd import DiscoveryClient
from nanover.trajectory import FRAME_SERVICE_NAME


def handle_user_arguments(args=None) -> argparse.Namespace:
    """
    Parse the arguments from the command line.

    :return: The namespace of arguments read from the command line.
    """
    description = textwrap.dedent(
        """\
    Recording client for NanoVer
    """
    )
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument(
        "--path",
        default=".",
        metavar="PATH",
        help="Destination directory for recording files.",
    )

    parser.add_argument(
        "--prefix",
        default="recording",
        metavar="PREFIX",
        help="Prefix for the recording filenames.",
    )

    parser.add_argument(
        "--address",
        metavar="ADDRESS:PORT",
        help="Connect to specific host address and port.",
    )

    parser.add_argument(
        "--autoconnect",
        metavar="QUERY",
        help="Use server discovery and pick a server with name matching the query.",
    )

    arguments = parser.parse_args(args)
    return arguments


def main():
    """
    Entry point for the command line.
    """

    arguments = handle_user_arguments()

    address = "localhost:38801"

    if arguments.address is not None:
        address = arguments.address
    elif arguments.autoconnect is not None:
        print(f"Using server discovery to find '{arguments.autoconnect}'")
        with DiscoveryClient() as discovery_client:
            for hub in discovery_client.search_for_services():
                if arguments.autoconnect in hub.name:
                    host, port = hub.get_service_address(FRAME_SERVICE_NAME)
                    address = f"{host}:{port}"
                    print(f"Found '{hub.name}' at {address}.")
                    break

    path = arguments.path
    prefix = arguments.prefix

    count = 0

    while True:
        count += 1
        timestamp = time.strftime("%Y%m%d-%H%M%S")
        name = f"{prefix}-{count}-{timestamp}"

        traj = f"{path}/{name}.traj"
        state = f"{path}/{name}.state"

        print(f"Connecting to {address}.")

        executor, channel = record_from_server(address, traj, state)

        print(f"Recording from server to {traj} and {state}. Press Ctrl-C to stop.")

        try:
            while True:
                time.sleep(0.01)
        except KeyboardInterrupt:
            pass
        finally:
            channel.close()
            executor.shutdown()

        try:
            input("Done. Press Ctrl-C to quit or Return to begin recording again.")
        except KeyboardInterrupt:
            break


if __name__ == "__main__":
    main()
