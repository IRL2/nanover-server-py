"""
Manually advertise a server on a web discovery endpoint.
"""

import argparse
import json
import textwrap

from nanover.app.discovery import DiscoveryEndpoint
from nanover.utilities.cli import suppress_keyboard_interrupt_as_cancellation


def handle_user_arguments(args=None) -> argparse.Namespace:
    """
    Parse the arguments from the command line.

    :return: The namespace of arguments read from the command line.
    """
    description = textwrap.dedent(
        """\
    Manually advertise a server on a web discovery endpoint.
    """
    )
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument(
        "--endpoint",
        metavar="HOST",
        help="Discovery endpoint to advertise on.",
    )

    parser.add_argument(
        "--name",
        default=None,
        metavar="NAME",
        help="Name of server to advertise.",
    )

    parser.add_argument(
        "--services",
        action="append",
        nargs="+",
        default=[],
        metavar="PROTOCOL://HOST[:PORT]",
        help="URLs of services to advertise.",
    )

    parser.add_argument(
        "--json",
        default="{}",
        metavar="JSON",
        help="Raw json data to advertise.",
    )

    arguments = parser.parse_args(args)
    return arguments


def main():
    """
    Entry point for the command line.
    """

    arguments = handle_user_arguments()
    endpoint = DiscoveryEndpoint.from_host(arguments.endpoint)

    data = json.loads(arguments.json)
    if arguments.name is not None:
        data["name"] = arguments.name
    for urls in arguments.services:
        for url in urls:
            data[url.split("://")[0]] = url

    with (
        suppress_keyboard_interrupt_as_cancellation() as cancellation,
        endpoint.advertise_data(data) as advert,
    ):
        print(advert)
        cancellation.wait_cancellation(interval=0.01)


if __name__ == "__main__":
    main()
