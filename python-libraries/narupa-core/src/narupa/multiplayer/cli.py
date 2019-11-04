# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.

import argparse
import textwrap
import time

from logging import StreamHandler

from narupa.essd import DiscoveryServer
from narupa.essd.servicehub import ServiceHub
from narupa.multiplayer.multiplayer_server import MultiplayerServer

"""
Command line interface for running a Narupa multiplayer server.
Run with:

.. code:: bash
    python cli.py

If the module is installed with pip, run with:

.. code:: bash
    narupa-multiplayer

"""


def handle_user_arguments() -> argparse.Namespace:
    """
    Parse the arguments from the command line.

    :return: The namespace of arguments read from the command line.
    """
    description = textwrap.dedent("""\
    Provide a multiplayer service with no other features.
    """)
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument(
        '-n', '--name',
        type=str, default='Narupa Multiplayer Server',
        help='Give a friendly name to the server.'
    )
    parser.add_argument('-p', '--port', default=54323)
    parser.add_argument('-a', '--address', default='[::]')
    parser.add_argument('-v', '--verbose', action="store_true", default=False)
    parser.add_argument('-vv', '--debug', action="store_true", default=False)
    parser.add_argument(
        '--no-discovery', dest='discovery', action='store_false', default=True,
        help='Run without the discovery service, so this server will not broadcast itself on the LAN.'
    )
    parser.add_argument(
        '--discovery-port', type=int, default=None,
        help='Port at which to run discovery service'
    )
    arguments = parser.parse_args()
    return arguments


def setup_discovery(server, name, discovery_port) -> DiscoveryServer:
    """
    Sets up an ESSD :class:`DiscoveryServer` for this multiplayer server.

    :param server: Multiplayer server.
    :param name: Name to broadcast with.
    :param discovery_port: Port on which to broadcast
    :return: DiscoveryServer with the multiplayer service registered and broadcasting.
    """
    discovery_server = DiscoveryServer(broadcast_port=discovery_port)
    service = ServiceHub(name=name, address=server.address)
    service.add_service(name="multiplayer", port=server.port)
    discovery_server.register_service(service)
    return discovery_server


def main():
    """
    Entry point for the command line.
    """
    arguments = handle_user_arguments()

    server = MultiplayerServer(address=arguments.address, port=arguments.port)

    if arguments.verbose:
        server._multiplayer_service.logger.setLevel("INFO")
    if arguments.debug:
        server._multiplayer_service.logger.setLevel("DEBUG")
    if arguments.verbose or arguments.debug:
        server._multiplayer_service.logger.addHandler(StreamHandler())

    print(f'Serving multiplayer on port {server.port}')

    discovery_server = None
    if arguments.discovery:
        discovery_server = setup_discovery(server, arguments.name, arguments.discovery_port)
    try:
        while True:
            time.sleep(1)
    except KeyboardInterrupt:
        print('Closing due to keyboard interrupt')
    finally:
        print('Cleaning up...')
        if discovery_server is not None:
            discovery_server.close()
        server.close()


if __name__ == '__main__':
    main()
