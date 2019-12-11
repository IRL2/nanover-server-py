# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.
"""
Interactive molecular dynamics runner for ASE with OpenMM.
"""
import logging
from typing import Optional

from ase import units, Atoms
from ase.md import MDLogger
from ase.md.nvtberendsen import NVTBerendsen
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from attr import dataclass
from narupa.ase import ase_to_frame_data, TrajectoryLogger
from narupa.ase.converter import add_ase_positions_to_frame_data
from narupa.ase.imd_server import ASEImdServer
from narupa.ase.wall_calculator import VelocityWallCalculator
from narupa.ase.openmm.calculator import OpenMMCalculator
from narupa.essd import DiscoveryServer
from narupa.essd.servicehub import ServiceHub
from narupa.imd.imd_service import IMD_SERVICE_NAME
from narupa.multiplayer import MultiplayerServer
from narupa.multiplayer.multiplayer_service import MULTIPLAYER_SERVICE_NAME
from narupa.openmm import openmm_to_frame_data, serializer
from narupa.core import get_requested_port_or_default
from narupa.trajectory.frame_publisher import TRAJECTORY_SERVICE_NAME
from narupa.trajectory.frame_server import DEFAULT_PORT as TRAJ_DEFAULT_PORT
from narupa.imd.imd_server import DEFAULT_PORT as IMD_DEFAULT_PORT
from narupa.multiplayer.multiplayer_server import DEFAULT_PORT as MULTIPLAYER_DEFAULT_PORT
from simtk.openmm.app import Simulation

CONSTRAINTS_UNSUPPORTED_MESSAGE = (
    "The simulation contains constraints which will be ignored by this runner!")


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
    name: str = None
    multiplayer: bool = True
    multiplayer_port: int = None
    discovery: bool = True
    discovery_port: int = None


@dataclass
class LoggingParams:
    """
    Class representing parameters for trajectory logging
    """
    trajectory_file: str = None
    log_interval: int = 1


class TrajectoryLoggerInfo:
    """
    Class giving a view into an ASE MD runners logger parameters.

    :param trajectory_logger: Trajectory logger performing the logging.
    :param params: Logging parameters.
    """

    def __init__(self, trajectory_logger: TrajectoryLogger, params:LoggingParams):
        """


        """
        self._logger = trajectory_logger
        self._params = params

    @property
    def trajectory_path(self) -> str:
        """
        The current trajectory path being outputted to.

        :return: The current trajectory path.
        """
        return self._logger.current_path

    @property
    def base_trajectory_path(self):
        """
        The base trajectory path, without timestamps.

        :return: The base trajectory path.
        """
        return self._logger.base_path

    @property
    def log_interval(self):
        """
        The interval at which logging is occurring.

        :return: The interval at which logging is occurring, in steps.
        """
        return self._params.log_interval


class OpenMMIMDRunner:
    """
    A wrapper class for running an interactive OpenMM simulation with ASE.

    :param simulation OpenMM simulation to run interactively.
    :param params IMD parameters to tune the server.
    :param logging_params Parameters for logging the trajectory of the simulation.
    """

    def __init__(self, simulation: Simulation,
                 imd_params: Optional[ImdParams] = None,
                 logging_params: Optional[LoggingParams] = None):
        self._logger = logging.getLogger(__name__)
        self.simulation = simulation
        self._validate_simulation()
        if not imd_params:
            imd_params = ImdParams()
        if not logging_params:
            logging_params = LoggingParams()
        self._address = imd_params.address
        if self._services_use_same_port(
                trajectory_port=imd_params.trajectory_port,
                imd_port=imd_params.imd_port,
                multiplayer_port=imd_params.multiplayer_port
        ):
            raise ValueError(
                "Trajectory serving port, IMD serving port and multiplayer serving port must be different!")
        self._frame_interval = imd_params.frame_interval
        self._time_step = imd_params.time_step
        self._verbose = imd_params.verbose

        self._initialise_calculator(simulation, walls=imd_params.walls)
        self._initialise_dynamics()
        self._initialise_server(self.dynamics,
                                imd_params.trajectory_port,
                                imd_params.imd_port,
                                imd_params.multiplayer,
                                imd_params.multiplayer_port,
                                imd_params.name,
                                imd_params.discovery,
                                imd_params.discovery_port)

        self._initialise_trajectory_logging(logging_params)

    def _validate_simulation(self):
        """
        Check this runner's simulation for unsupported features and issue the
        relevant warnings.
        """
        if self.simulation.system.getNumConstraints() > 0:
            self._logger.warning(CONSTRAINTS_UNSUPPORTED_MESSAGE)

    @classmethod
    def from_xml(cls, simulation_xml,
                 params: Optional[ImdParams] = None,
                 logging_params: Optional[LoggingParams] = None):
        """
        Initialises a :class:`OpenMMIMDRunner` from a simulation XML file
        serialised with :fun:`serializer.serialize_simulation`.

        :param simulation_xml: Path to XML file.
        :param params: The :class: ImdParams to run the server with.
        :param logging_params: The :class:LoggingParams to set up trajectory logging with.
        :return: An OpenMM simulation runner.
        """
        with open(simulation_xml) as infile:
            simulation = serializer.deserialize_simulation(infile.read())
        return OpenMMIMDRunner(simulation, params, logging_params)

    @property
    def verbose(self):
        """
        Whether this OpenMM runner is set to run in verbose mode. If it is, it
        will print state information every 100 steps using an :class: MDLogger.

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
        Gets the port the :class:`FrameServer` is running on.

        :return: The port the frame service is running on.
        """
        return self.imd.frame_server.port

    @property
    def imd_port(self):
        """
        Gets the port the :class:`ImdServer` is running on.

        :return: The port the IMD server is running on.
        """
        return self.imd.imd_server.port

    @property
    def running_multiplayer(self):
        return self.multiplayer is not None

    @property
    def multiplayer_port(self):
        try:
            return self.multiplayer.port
        except AttributeError:
            raise AttributeError("Multiplayer service not running")

    @property
    def name(self):
        return self.imd.name

    @property
    def running_discovery(self):
        try:
            return self.discovery_server is not None
        except AttributeError:
            return False

    @property
    def discovery_port(self):
        try:
            return self.discovery_server.port
        except AttributeError:
            raise AttributeError("Discovery service not running")

    @property
    def dynamics(self):
        """
        Gets the ASE :class:`MolecularDynamics` object that is running the dynamics.

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
        if self.multiplayer is not None:
            self.multiplayer.close()
        if self.discovery_server is not None:
            self.discovery_server.close()

    def _initialise_server(self, dynamics,
                           trajectory_port=None,
                           imd_port=None,
                           run_multiplayer=True,
                           multiplayer_port=None,
                           name=None,
                           run_discovery=True,
                           discovery_port=None):
        # set the server to use the OpenMM frame convert for performance purposes.
        self.imd = ASEImdServer(dynamics,
                                frame_method=openmm_ase_frame_server,
                                address=self.address,
                                frame_interval=self.frame_interval,
                                trajectory_port=trajectory_port,
                                imd_port=imd_port,
                                name=name,
                                )

        if run_multiplayer:
            self.multiplayer = MultiplayerServer(address=self.address, port=multiplayer_port)
        else:
            self.multiplayer = None

        if run_discovery:
            self.discovery_server = DiscoveryServer(discovery_port)
            self._register_services(name)
        else:
            self.discovery_server = None

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

    def _initialise_trajectory_logging(self, logging_params: LoggingParams):
        if not logging_params.trajectory_file:
            self.logging_info = None
            return
        logger = TrajectoryLogger(self.atoms, logging_params.trajectory_file)
        self.dynamics.attach(logger, logging_params.log_interval)
        self.imd.on_reset_listeners.append(logger.reset)
        self.logging_info = TrajectoryLoggerInfo(logger, logging_params)

    @staticmethod
    def _services_use_same_port(trajectory_port, imd_port, multiplayer_port):
        trajectory_port = get_requested_port_or_default(trajectory_port, TRAJ_DEFAULT_PORT)
        imd_port = get_requested_port_or_default(imd_port, IMD_DEFAULT_PORT)
        multiplayer_port = get_requested_port_or_default(multiplayer_port, MULTIPLAYER_DEFAULT_PORT)
        # If a port is set to 0, then GRPC will choose one available port; so
        # 0 is always a valid value.
        all_ports = (trajectory_port, imd_port, multiplayer_port)
        non_zero_ports = [port for port in all_ports if port != 0]
        return len(non_zero_ports) != len(set(non_zero_ports))

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        self.close()

    def _register_services(self, server_name):
        hub = ServiceHub(name=server_name, address=self.imd.frame_server.address)
        hub.add_service(name=IMD_SERVICE_NAME, port=self.imd.imd_server.port)
        hub.add_service(name=TRAJECTORY_SERVICE_NAME, port=self.imd.frame_server.port)
        if self.multiplayer is not None:
            hub.add_service(name=MULTIPLAYER_SERVICE_NAME, port=self.multiplayer.port)
        self.discovery_server.register_service(hub)
