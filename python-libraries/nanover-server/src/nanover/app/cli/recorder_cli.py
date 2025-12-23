"""
A NanoVer client that records the trajectory and the state updates.

Run the script as `nanover-record` to connect to a server running on localhost
on the default port (38801) and begin recording to uniquely named .traj and
.state files. Run `nanover-record --help` to find options for connecting to
specific servers and changing the output filename.
"""

import argparse
import textwrap
import time

from nanover.essd import DiscoveryClient
from nanover.utilities.cli import suppress_keyboard_interrupt_as_cancellation
from nanover.websocket.client.app_client import get_websocket_address_from_hub
from nanover.websocket.record import BackgroundRecordingContext


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
        help="Destination for recording files.",
    )

    parser.add_argument(
        "--prefix",
        default="recording",
        metavar="PREFIX",
        help="Prefix for the recording filenames.",
    )

    parser.add_argument(
        "--address",
        metavar="PROTOCOL://HOST:PORT",
        help="Connect to specific host address and port.",
        default="ws://127.0.0.1:38801",
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
    address = arguments.address

    if arguments.autoconnect is not None:
        print(f"Using server discovery to find '{arguments.autoconnect}'")
        with DiscoveryClient() as discovery_client:
            for hub in discovery_client.search_for_services():
                if arguments.autoconnect in hub.name:
                    address = get_websocket_address_from_hub(hub)
                    print(f"Found '{hub.name}' at {address}.")
                    break

    path = arguments.path
    prefix = arguments.prefix

    count = 0

    while True:
        count += 1
        timestamp = time.strftime("%Y%m%d-%H%M%S")
        name = f"{prefix}-{count}-{timestamp}"

        outfile = f"{path}/{name}.nanover.zip"

        print(f"Connecting to {address}.")

        interrupted = False

        with suppress_keyboard_interrupt_as_cancellation() as cancellation:
            with BackgroundRecordingContext.from_address_to_path(
                address=address, path=outfile
            ) as recording:
                print(f"Recording from server to {outfile}. Press Ctrl-C to stop.")
                recording.future.add_done_callback(lambda _: cancellation.cancel())
                cancellation.wait_cancellation(interval=0.01)
                interrupted = not recording.future.done()

        if not interrupted:
            print("Recording ended by disconnection.")
            break

        try:
            input(
                "Recording ended by Ctrl-C. Press Ctrl-C to quit or Return to begin recording again."
            )
        except KeyboardInterrupt:
            break
        except EOFError:
            break


if __name__ == "__main__":
    main()
