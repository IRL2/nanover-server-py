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
    parser.add_argument('-t', '--trajectory_port', type=int, default=None)
    parser.add_argument('-i', '--imd_port', type=int, default=None)
    parser.add_argument('-a', '--address', default=None)
    parser.add_argument('-f', '--frame_interval', type=int, default=5)
    parser.add_argument('-s', '--time_step', type=float, default=2.0)
    arguments = parser.parse_args(args)
    return arguments


def initialise(args=None):
    arguments = handle_user_arguments(args)

    # TODO clean way to handle params?
    params = ImdParams(arguments.address,
                       arguments.trajectory_port,
                       arguments.imd_port,
                       arguments.frame_interval,
                       arguments.time_step,
                       arguments.verbose)
    runner = OpenMMIMDRunner.from_xml(arguments.simulation_xml_path, params)
    return runner


def main():
    """
    Entry point for the command line.
    """
    with initialise() as runner:
        print(f'Serving frames on port {runner.trajectory_port} and IMD on {runner.imd_port}')
        
        try:
            while True:
                runner.run(100)
        except KeyboardInterrupt:
            print("Closing due to keyboard interrupt.")

if __name__ == '__main__':
    main()
