# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.
"""
Interactive molecular dynamics runner for ASE with OpenMM.
"""
from typing import Optional

from ase import units, Atoms
from ase.md import MDLogger
from ase.md.nvtberendsen import NVTBerendsen
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from attr import dataclass
from narupa.ase import ase_to_frame_data
from narupa.ase.converter import add_ase_positions_to_frame_data
from narupa.ase.imd_server import ASEImdServer
from narupa.ase.openmm.calculator import OpenMMCalculator
from narupa.openmm import openmm_to_frame_data, serializer
from simtk.openmm.app import Simulation


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


@dataclass
class ImdParams:
    """
    Class representing parameters for IMD runners.
    """
    address: str = None
    trajectory_port: int = None
    imd_port: int = None
    frame_interval: int = 5
    time_step: float = 1.0
    verbose: bool = False


class OpenMMIMDRunner:
    """
    A wrapper class for running an interactive OpenMM simulation with ASE.

    :param simulation OpenMM simulation to run interactively.
    :param params IMD parameters to tune the server.
    """
    def __init__(self, simulation:Simulation, params: Optional[ImdParams] = None):
        self.simulation = simulation
        if not params:
            params = ImdParams()
        self._address = params.address
        self._trajectory_port = params.trajectory_port
        self._imd_port = params.imd_port
        self._frame_interval = params.frame_interval
        self._time_step = params.time_step
        self._verbose = params.verbose

        self._initialise_calculator(simulation)
        self._initialise_dynamics(simulation)
        self._initialise_server(self.dyn)

    @classmethod
    def from_xml(cls, simulation_xml, params: Optional[ImdParams] = None):
        """
        Initialises a :class: OpenMMIMDRunner from a simulation XML file serialised with
        :class: serializer
        :param simulation_xml: Path to XML file.
        :param params: The :class: ImdParams to run the server with.
        :return: An OpenMM simulation runner.
        """
        with open(simulation_xml) as infile:
            simulation = serializer.deserialize_simulation(infile.read())
        return OpenMMIMDRunner(simulation, params)

    @property
    def verbose(self):
        return self._verbose

    @property
    def time_step(self):
        return self._time_step

    @property
    def frame_interval(self):
        return self._frame_interval

    @property
    def address(self):
        return self._address

    @property
    def trajectory_port(self):
        return self._trajectory_port

    @property
    def imd_port(self):
        return self._imd_port

    def run(self, steps=None):
        self.imd.run(steps)

    def _initialise_server(self, dynamics):
        # set the server to use the OpenMM frame convert for performance purposes.
        self.imd = ASEImdServer(dynamics,
                                frame_method=openmm_ase_frame_server,
                                address=self.address,
                                frame_interval=self.frame_interval,
                                trajectory_port=self.trajectory_port,
                                imd_port=self.imd_port,
                                )

    def _initialise_calculator(self, simulation):
        self._openmm_calculator = OpenMMCalculator(simulation)
        self.atoms = self._openmm_calculator.generate_atoms()
        self.atoms.set_calculator(self._openmm_calculator)

    def _initialise_dynamics(self, simulation):
        # Set the momenta corresponding to T=300K
        MaxwellBoltzmannDistribution(self.atoms, 300 * units.kB)

        self.dyn = NVTBerendsen(self.atoms, 1 * units.fs, 300, self.time_step * units.fs)

        if self.verbose:
            self.dyn.attach(MDLogger(self.dyn, self.atoms, '-', header=True, stress=False,
                                     peratom=False), interval=100)
