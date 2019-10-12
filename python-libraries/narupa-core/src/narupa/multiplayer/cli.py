# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.

import argparse
import textwrap
import time

from logging import StreamHandler

from narupa.multiplayer.multiplayer_server import MultiplayerServer

"""
Command line interface for running a Narupa multiplayer server.
Run with:

.. code bash
    python cli.py

If the module is installed with pip, run with:
.. code bash
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

    parser.add_argument('-p', '--port', default=54323)
    parser.add_argument('-a', '--address', default='[::]')
    parser.add_argument('-v', '--verbose', action="store_true", default=False)
    parser.add_argument('-vv', '--debug', action="store_true", default=False)
    arguments = parser.parse_args()
    return arguments


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

    try:
        while True:
            time.sleep(1)
    except KeyboardInterrupt:
        print('Closing due to keyboard interrupt')
    finally:
        server.close()


if __name__ == '__main__':
    main()
