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

from narupa.app.app_server import qualified_server_name
from narupa.ase.openmm import ASEOpenMMRunner
from narupa.ase.openmm.runner import ImdParams, LoggingParams


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
        type=str, default=None,
        help='Give a friendly name to the server.'
    )
    parser.add_argument('-p', '--port', type=int, default=None)
    parser.add_argument('-a', '--address', default=None)
    parser.add_argument(
        '-f', '--frame-interval', type=int, default=5,
        help='Produce a trajectory frame every LOG_INTERVAL dynamics steps.'
    )
    parser.add_argument(
        '-s', '--time-step', type=float, default=1.0,
        help='The simulation time step, in femtoseconds.',
    )
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
    parser.add_argument(
        '-o', '--trajectory-file', type=str, default=None,
        help='Base filename for the trajectory output file. A timestamp will '
             'be inserted between the name and file extensions. Can be any '
             'file that ASE can output in append mode, such as XYZ.'
    )
    parser.add_argument(
        '--write-interval', type=int, default=1,
        help='Write a trajectory frame to file every WRITE_INTERVAL dynamics '
             'steps.',
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

    if arguments.name is None:
        arguments.name = qualified_server_name("Narupa OpenMM ASE Server")

    logging_params = LoggingParams(
        arguments.trajectory_file,
        arguments.write_interval,
    )
    runner = ASEOpenMMRunner.from_xml(arguments.simulation_xml_path, params, logging_params)
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
        if runner.logging_info is not None:
            print(f'Logging frames to "{runner.logging_info.trajectory_path}"')

        runner.imd.on_reset_listeners.append(lambda: print('RESET! ' * 10))
        print(f'Serving "{runner.name}" on port {runner.app_server.port}, discoverable on all interfaces on port {runner.discovery_port}')
        try:
            runner.run(block=False, reset_energy=runner.cli_options['reset_energy'])
            while True:
                time.sleep(1)
        except KeyboardInterrupt:
            print("Closing due to keyboard interrupt.")


if __name__ == '__main__':
    main()
