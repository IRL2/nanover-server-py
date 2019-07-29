"""
Command line interface for connecting to a narupa server and OSC server and
generating outgoing OSC messages from incoming narupa frame data. Used f
"""
import argparse
import textwrap

from osc_client import OscClient


class OscApp:
    def __init__(self, setup_callback=None):
        self._argument_parser = self._create_argument_parser()
        self._setup_callback = setup_callback

    def _create_argument_parser(self):
        description = textwrap.dedent("""\
                    Connect to a narupa trajectory service, and osc server and generate outgoing
                    osc messages from incoming frame data.
                    """)
        parser = argparse.ArgumentParser(description=description)

        parser.add_argument('--traj-address', default=None)
        parser.add_argument('--osc-address', default='localhost')
        parser.add_argument('-t', '--traj-port', type=int, default=None)
        parser.add_argument('-o', '--osc-port', type=int, default=None)
        parser.add_argument('-i', '--send-interval', type=float, default=.01)
        parser.add_argument('-v', '--verbose', action="store_true", default=False)
        return parser

    def _create_client(self, args):
        arguments = self._argument_parser.parse_args(args)
        osc_client = OscClient(traj_address=arguments.traj_address,
                               traj_port=arguments.traj_port,
                               osc_address=arguments.osc_address,
                               osc_port=arguments.osc_port,
                               send_interval=arguments.send_interval,
                               verbose=arguments.verbose)
        return osc_client

    def run(self, args=None):
        with self._create_client(args) as osc_client:
            self._setup_callback(osc_client)
            print('Running...')
            try:
                osc_client.run()
            except KeyboardInterrupt:
                print("Closing due to keyboard interrupt.")
