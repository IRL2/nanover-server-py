"""
Command line interface for running an interactive OpenMM simulation with ASE.
Run with:

.. code bash
    python cli.py neuraminidase.xml

If the module is installed with pip, run with:

.. code bash
    nanover-omm-ase neuraminidase.xml

"""

import argparse
import textwrap
from contextlib import contextmanager

from nanover.app import NanoverImdApplication
from nanover.omni import OmniRunner
from nanover.omni.ase_omm import ASEOpenMMSimulation


def handle_user_arguments(args=None) -> argparse.Namespace:
    """
    Parse the arguments from the command line.

    :return: The namespace of arguments read from the command line.
    """
    description = textwrap.dedent(
        """\
    Run an ASE IMD simulation with an OpenMM simulation.
    """
    )
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument(
        "simulation_xml_paths",
        nargs="+",
        help="The simulations to run in XML format.",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Display state information.",
    )
    parser.add_argument(
        "-n",
        "--name",
        type=str,
        default=None,
        help="Give a friendly name to the server.",
    )
    parser.add_argument("-p", "--port", type=int, default=None)
    parser.add_argument("-a", "--address", default=None)
    parser.add_argument(
        "-f",
        "--frame-interval",
        type=int,
        default=5,
        help="Produce a trajectory frame every LOG_INTERVAL dynamics steps.",
    )
    parser.add_argument(
        "-s",
        "--time-step",
        type=float,
        default=1.0,
        help="The simulation time step, in femtoseconds.",
    )
    parser.add_argument(
        "--reset-energy",
        type=float,
        default=1e6,
        help=(
            "Threshold of total energy above which the simulation is reset "
            "(kJ/mol). The value is ignored if --no-auto-reset is used."
        ),
    )
    parser.add_argument(
        "--no-auto-reset",
        dest="auto_reset",
        action="store_false",
        default=True,
        help="Do not reset the simulation, even if the energy becomes high.",
    )
    parser.add_argument(
        "-w",
        "--walls",
        action="store_true",
        default=False,
        help="Set a wall around the box, atoms will bounce against it.",
    )
    parser.add_argument(
        "--platform", default=None, help="Select the platform on which to run Openmm."
    )
    arguments = parser.parse_args(args)
    return arguments


@contextmanager
def initialise(args=None):
    arguments = handle_user_arguments(args)

    app_server = NanoverImdApplication.basic_server(
        name=arguments.name,
        address=arguments.address,
        port=arguments.port,
    )
    runner = OmniRunner(app_server)

    for path in arguments.simulation_xml_paths:
        simulation = ASEOpenMMSimulation.from_xml_path(path)
        simulation.verbose = arguments.verbose
        simulation.platform = arguments.platform
        simulation.frame_interval = arguments.frame_interval
        simulation.time_step = arguments.time_step
        simulation.use_walls = arguments.walls
        simulation.reset_energy = (
            arguments.reset_energy if arguments.auto_reset else None
        )

        simulation.on_reset_energy_exceeded.add_callback(lambda: print("RESET! " * 10))

        runner.add_simulation(simulation)

    yield runner


def main():
    """
    Entry point for the command line.
    """

    with initialise() as runner:
        runner.print_basic_info_and_wait()


if __name__ == "__main__":
    main()
