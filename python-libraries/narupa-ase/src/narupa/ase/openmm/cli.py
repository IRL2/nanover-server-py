"""
Command line interface for running an interactive OpenMM simulation with ASE.
Run with:

.. code bash
    python cli.py neuraminidase.xml

If the module is installed with pip, run with:
.. code bash
    narupa-omm-ase neuraminidase.xml

"""
import argparse
import textwrap
import time

from narupa.ase.openmm import OpenMMIMDRunner
from narupa.ase.openmm.runner import ImdParams


def handle_user_arguments(args=None) -> argparse.Namespace:
    """
    Parse the arguments from the command line.

    :return: The namespace of arguments read from the command line.
    """
    description = textwrap.dedent("""\
    Run an ASE IMD simulation with an OpenMM simulation.
    """)
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument(
        'simulation_xml_path',
        help='The simulation to run in XML format.',
    )
    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='Display state information.',
    )
    parser.add_argument(
        '-n', '--name',
        type=str, default='Narupa OpenMM ASE Server',
        help='Give a friendly name to the server.'
    )
    parser.add_argument('-p', '--port', type=int, default=None)
    parser.add_argument('-a', '--address', default=None)
    parser.add_argument('-f', '--frame-interval', type=int, default=5)
    parser.add_argument('-s', '--time-step', type=float, default=1.0)
    parser.add_argument(
        '--reset-energy', type=float, default=1e6,
        help=('Threshold of total energy above which the simulation is reset '
              '(kJ/mol). The value is ignored if --no-auto-reset is used.'),
    )
    parser.add_argument(
        '--no-auto-reset', dest='auto_reset', action='store_false', default=True,
        help='Do not reset the simulation, even if the energy becomes high.',
    )
    parser.add_argument(
        '-w', '--walls', action='store_true', default=False,
        help='Set a wall around the box, atoms will bounce against it.',
    )
    parser.add_argument(
        '--no-discovery', dest='discovery', action='store_false', default=True,
        help='Run without the discovery service, so this server will not broadcast itself on the LAN.'
    )
    parser.add_argument(
        '--discovery-port', type=int, default=None,
        help='Port at which to run discovery service'
    )
    arguments = parser.parse_args(args)
    return arguments


def initialise(args=None):
    arguments = handle_user_arguments(args)

    # TODO clean way to handle params?
    params = ImdParams(
        arguments.address,
        arguments.port,
        arguments.frame_interval,
        arguments.time_step,
        arguments.verbose,
        arguments.walls,
        arguments.name,
        arguments.discovery,
        arguments.discovery_port
    )
    runner = OpenMMIMDRunner.from_xml(arguments.simulation_xml_path, params)
    # Shamefully store CLI arguments in the runner.
    runner.cli_options = {
        'reset_energy': arguments.reset_energy if arguments.auto_reset else None,
    }
    return runner


def main():
    """
    Entry point for the command line.
    """
    with initialise() as runner:
        runner.imd.on_reset_listeners.append(lambda: print('RESET! ' * 10))
        print(f'Serving on port {runner.server.port}, discoverable on all interfaces on port {runner.discovery_port}')
        try:
            runner.run(block=False, reset_energy=runner.cli_options['reset_energy'])
            while True:
                time.sleep(1)
        except KeyboardInterrupt:
            print("Closing due to keyboard interrupt.")


if __name__ == '__main__':
    main()
