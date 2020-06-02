# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.
import argparse
import textwrap

from narupa.app import NarupaImdClient
from osc_client import OscClient, DEFAULT_OSC_ADDRESS


class OscApp:
    """
    For connecting to a Narupa simulation and then uses the result of the
    message generator setup to translate incoming trajectory frames into
    outgoing osc messages.
    """
    def __init__(self, message_generator_setup):
        self._argument_parser = self._create_argument_parser()
        self._message_generator_setup = message_generator_setup

    @staticmethod
    def _create_argument_parser():
        description = textwrap.dedent("""\
                    Connect to a narupa trajectory service, and osc server and generate outgoing
                    osc messages from incoming frame data.
                    """)
        parser = argparse.ArgumentParser(description=description)

        default_host, default_port = DEFAULT_OSC_ADDRESS
        parser.add_argument('--traj-address', default=None)
        parser.add_argument('--osc-host', default=default_host)
        parser.add_argument('-n', '--server-name', type=str, default=None)
        parser.add_argument('-o', '--osc-port', type=int, default=default_port)
        parser.add_argument('-i', '--send-interval', type=float, default=.01)
        parser.add_argument('-v', '--verbose', action="store_true", default=False)
        return parser

    def _create_client(self, args):
        arguments = self._argument_parser.parse_args(args)
        osc_client = OscClient(narupa_client=NarupaImdClient.autoconnect(),
                               osc_address=(arguments.osc_host, arguments.osc_port),
                               osc_send_interval=arguments.send_interval,
                               verbose=arguments.verbose)
        osc_client.narupa_client.subscribe_multiplayer()
        return osc_client

    def run(self, args=None):
        """
        Begin connection, setup, and then run until keyboard interrupt.
        """
        server_name = args.server_name if args is not None else None

        with self._create_client(args) as osc_client:
            osc_client.message_generator = self._message_generator_setup(osc_client)
            print('Running...')
            try:
                osc_client.run()
            except KeyboardInterrupt:
                print("Closing due to keyboard interrupt.")
