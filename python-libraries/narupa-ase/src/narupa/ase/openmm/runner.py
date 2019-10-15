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
from narupa.ase.wall_calculator import VelocityWallCalculator
from narupa.ase.openmm.calculator import OpenMMCalculator
from narupa.openmm import openmm_to_frame_data, serializer
from narupa.core import get_requested_port_or_default
from narupa.trajectory.frame_server import DEFAULT_PORT as TRAJ_DEFAULT_PORT
from narupa.imd.imd_server import DEFAULT_PORT as IMD_DEFAULT_PORT
from simtk.openmm.app import Simulation


def openmm_ase_frame_server(ase_atoms: Atoms, frame_server):
    """
    Generates and sends frames for a simulation using an :class: OpenMMCalculator.
    """

    def send():
        # generate topology frame using OpenMM converter.
        if send.frame_index == 0:
            imd_calculator = ase_atoms.get_calculator()
            send.topology = imd_calculator.calculator.topology
            frame = openmm_to_frame_data(state=None, topology=send.topology)
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
    walls: bool = False


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
        if self._services_use_same_port(
                trajectory_port=params.trajectory_port,
                imd_port=params.imd_port,
            ):
            raise ValueError("Trajectory serving port and IMD serving port must be different!")
        self._frame_interval = params.frame_interval
        self._time_step = params.time_step
        self._verbose = params.verbose

        self._initialise_calculator(simulation, walls=params.walls)
        self._initialise_dynamics()
        self._initialise_server(self.dynamics, params.trajectory_port, params.imd_port)

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
        """
        Whether this OpenMM runner is set to run in verbose mode. If it is, it will print state information
        every 100 steps using an :class: MDLogger.
        :return: `True` if set to run verbosely, `False` otherwise.
        """
        return self._verbose

    @property
    def time_step(self):
        """
        Gets the time step of the simulation, in femtoseconds.
        :return: The time step of the simulation.
        """
        return self._time_step

    @property
    def frame_interval(self):
        """
        Gets the interval at which frames are sent, in steps.
        :return: The frame interval, in steps.
        """
        return self._frame_interval

    @property
    def address(self):
        """
        Gets the URL or IP address the server is running at.
        :return: The URL or IP address of the server.
        """
        return self._address

    @property
    def trajectory_port(self):
        """
        Gets the port the :class: FrameServer is running on.
        :return: The port the frame service is running on.
        """
        return self.imd.frame_server.port

    @property
    def imd_port(self):
        """
        Gets the port the :class: ImdServer is running on.
        :return: The port the IMD server is running on.
        """
        return self.imd.imd_server.port

    @property
    def dynamics(self):
        """
        Gets the ASE :class: MolecularDynamics object that is running the dynamics.
        :return: The ASE molecular dynamics object.
        """
        return self._dynamics

    def run(self, steps: Optional[int] = None,
            block: Optional[bool] = None, reset_energy: Optional[float] = None):
        """
        Runs the molecular dynamics.

        :param steps: If passed, will run the given number of steps, otherwise
            will run forever on a background thread and immediately return.
        :param block: If ``False`` run in a separate thread. By default, "block"
            is ``None``, which means it is automatically set to ``True`` if a
            number of steps is provided and to ``False`` otherwise.
        :param reset_energy: Threshold of total energy in kJ/mol above which
            the simulation is reset to its initial conditions. If a value is
            provided, the simulation is reset if the total energy is greater
            than this value, or if the total energy is `nan` or infinite. If
            ``None`` is provided instead, then the simulation will not be
            automatically reset.
        """
        self.imd.run(steps, block=block, reset_energy=reset_energy)

    def close(self):
        """
        Closes the connection and stops the dynamics.
        """
        self.imd.close()

    def _initialise_server(self, dynamics, trajectory_port=None, imd_port=None):
        # set the server to use the OpenMM frame convert for performance purposes.
        self.imd = ASEImdServer(dynamics,
                                frame_method=openmm_ase_frame_server,
                                address=self.address,
                                frame_interval=self.frame_interval,
                                trajectory_port=trajectory_port,
                                imd_port=imd_port,
                                )

    def _initialise_calculator(self, simulation, walls=False):
        self._openmm_calculator = OpenMMCalculator(simulation)
        self._md_calculator = self._openmm_calculator
        self.atoms = self._openmm_calculator.generate_atoms()
        if walls:
            self._md_calculator = VelocityWallCalculator(self._openmm_calculator, self.atoms)
        self.atoms.set_calculator(self._md_calculator)

    def _initialise_dynamics(self):
        # Set the momenta corresponding to T=300K
        MaxwellBoltzmannDistribution(self.atoms, 300 * units.kB)

        # We do not remove the center of mass (fixcm=False). If the center of
        # mass translations should be removed, then the removal should be added
        # to the OpenMM system.
        self._dynamics = NVTBerendsen(
            atoms=self.atoms,
            timestep=self.time_step * units.fs,
            temperature=300,
            taut=self.time_step * units.fs,
            fixcm=False,
        )

        if self.verbose:
            self._dynamics.attach(MDLogger(self._dynamics, self.atoms, '-', header=True, stress=False,
                                     peratom=False), interval=100)
    
    @staticmethod
    def _services_use_same_port(trajectory_port, imd_port):
        trajectory_port = get_requested_port_or_default(trajectory_port, TRAJ_DEFAULT_PORT)
        imd_port = get_requested_port_or_default(imd_port, IMD_DEFAULT_PORT)
        # If a port is set to 0, then GRPC will affect one available port; so
        # 0 is always a valid value.
        if trajectory_port == 0 or imd_port == 0:
            return False
        return trajectory_port == imd_port

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        self.close()
