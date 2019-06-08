# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.
"""
Demonstrates interactive molecular dynamics running with OpenMM as the engine and ASE as the integrator.

Run with:

.. code bash
    python imd_openmm.py neuraminidase.xml

"""
import argparse
import textwrap

from ase import units, Atoms
from ase.md import MDLogger
from ase.md.nvtberendsen import NVTBerendsen
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from narupa.openmm import openmm_to_frame_data

from narupa.ase import ase_to_frame_data
from narupa.ase.converter import add_ase_positions_to_frame_data
from narupa.ase.imd_server import ASEImdServer
from openmm_calculator import OpenMMCalculator


def openmm_ase_frame_server(ase_atoms: Atoms, frame_server):
    """
    Generates and sends frames for a simulation using an :class: OpenMMCalculator.
    """

    def send():
        # generate topology frame using OpenMM converter.
        if send.frame_index == 0:
            imd_calculator = ase_atoms.get_calculator()
            send.topology = imd_calculator.calculator.simulation.topology
            frame = openmm_to_frame_data(positions=None, topology=send.topology)
            add_ase_positions_to_frame_data(frame, ase_atoms.get_positions())
        # from then on, just send positions and state.
        else:
            frame = ase_to_frame_data(ase_atoms, topology=False)
        frame_server.send_frame(send.frame_index, frame)
        send.frame_index = send.frame_index + 1

    send.frame_index = 0
    send.topology = None
    return send


def handle_user_arguments() -> argparse.Namespace:
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
    parser.add_argument('-t', '--trajectory_port', default=54321)
    parser.add_argument('-i', '--imd_port', default=54322)
    parser.add_argument('-a', '--address', default='[::]')
    parser.add_argument('-f', '--frame_interval', default=5)
    parser.add_argument('--time_step', default=2.0)
    arguments = parser.parse_args()
    return arguments


def main():
    """
    Entry point for the command line.
    """
    arguments = handle_user_arguments()

    input_xml = arguments.simulation_xml_path
    print(f'Generating OpenMM context from input: {input_xml}')
    openmm_calculator = OpenMMCalculator(input_xml)

    print(f'Generating ASE representation of the OpenMM system')
    atoms = openmm_calculator.generate_atoms()

    atoms.set_calculator(openmm_calculator)
    print(f'ASE energy: {atoms.get_potential_energy()}')
    # Set the momenta corresponding to T=300K
    MaxwellBoltzmannDistribution(atoms, 300 * units.kB)

    dyn = NVTBerendsen(atoms, 1 * units.fs, 300, arguments.time_step * units.fs)

    if arguments.verbose:
        dyn.attach(MDLogger(dyn, atoms, '-', header=True, stress=False,
                            peratom=False), interval=100)
    # set the server to use the OpenMM frame convert for performance purposes.
    imd = ASEImdServer(dyn,
                       frame_method=openmm_ase_frame_server,
                       address=arguments.address,
                       frame_interval=arguments.frame_interval,
                       trajectory_port=arguments.trajectory_port,
                       imd_port=arguments.imd_port,
                       )
    print(f'Running dynamics')
    while True:
        imd.run(100)


if __name__ == '__main__':
    main()
