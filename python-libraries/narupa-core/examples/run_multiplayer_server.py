# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.
"""
Demonstrates multiplayer server with no additional features.

Run with:

.. code bash
    python multiplayer_only.py

"""
import argparse
import textwrap
import time

from logging import StreamHandler

from narupa.multiplayer.multiplayer_server import MultiplayerServer

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
    parser.add_argument('-s', '--send-self', action="store_true", default=False)
    parser.add_argument('-v', '--verbose', action="store_true", default=False)
    parser.add_argument('-vv', '--debug', action="store_true", default=False)
    arguments = parser.parse_args()
    return arguments


def main():
    """
    Entry point for the command line.
    """
    arguments = handle_user_arguments()
    
    server = MultiplayerServer(address=arguments.address, port=arguments.port, send_self=arguments.send_self)
    
    if arguments.verbose:
        server.multiplayer_services.logger.setLevel("INFO")
        server.multiplayer_services.logger.addHandler(StreamHandler())
    if argument.debug:
        server.multiplayer_services.logger.setLevel("DEBUG")
        server.multiplayer_services.logger.addHandler(StreamHandler())

    print(f'Serving multiplayer on port {arguments.port}')

    try:
        while True:
            time.sleep(1)
    except KeyboardInterrupt:
        print('Good bye.')
    finally:
        server.close()

if __name__ == '__main__':
    main()
