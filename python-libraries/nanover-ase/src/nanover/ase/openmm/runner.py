"""
Interactive molecular dynamics runner for ASE with OpenMM.
"""

import logging
from typing import Optional, List

from ase import units, Atoms  # type: ignore
from ase.md import MDLogger, Langevin
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from attr import dataclass
from nanover.app import NanoverImdApplication, NanoverRunner
from nanover.app.app_server import DEFAULT_NANOVER_PORT
from nanover.core import NanoverServer, DEFAULT_SERVE_ADDRESS
from nanover.ase import TrajectoryLogger
from nanover.essd import DiscoveryServer
from nanover.openmm import openmm_to_frame_data, serializer
from nanover.trajectory.frame_publisher import FramePublisher
from simtk.openmm.app import Simulation, Topology

from nanover.ase import ase_to_frame_data
from nanover.ase.converter import add_ase_positions_to_frame_data
from nanover.ase.imd import NanoverASEDynamics
from nanover.ase.openmm.calculator import OpenMMCalculator
from nanover.ase.wall_constraint import VelocityWallConstraint

CONSTRAINTS_UNSUPPORTED_MESSAGE = (
    "The simulation contains constraints which will be ignored by this runner!"
)


def openmm_ase_frame_adaptor(ase_atoms: Atoms, frame_publisher: FramePublisher):
    """
    Generates and sends frames for a simulation using an :class: OpenMMCalculator.
    """

    frame_index = 0
    topology: Optional[Topology] = None

    def send():
        nonlocal frame_index, topology
        # generate topology frame using OpenMM converter.
        if frame_index == 0:
            imd_calculator = ase_atoms.calc
            topology = imd_calculator.calculator.topology
            frame = openmm_to_frame_data(state=None, topology=topology)
            add_ase_positions_to_frame_data(frame, ase_atoms.get_positions())
        # from then on, just send positions and state.
        else:
            frame = ase_to_frame_data(ase_atoms, topology=False)
        frame_publisher.send_frame(frame_index, frame)
        frame_index += 1

    return send


@dataclass
class ImdParams:
    """
    Class representing parameters for IMD runners.
    """

    address: Optional[str] = None
    port: Optional[int] = None
    frame_interval: int = 5
    time_step: float = 1.0
    verbose: bool = False
    walls: bool = False
    name: Optional[str] = None
    discovery: bool = True
    discovery_port: Optional[int] = None
    platform: Optional[str] = None


@dataclass
class LoggingParams:
    """
    Class representing parameters for trajectory logging
    """

    trajectory_file: Optional[str] = None
    write_interval: int = 1


class TrajectoryLoggerInfo:
    """
    Class giving a view into an ASE MD runners logger parameters.

    :param trajectory_logger: Trajectory logger performing the logging.
    :param params: Logging parameters.
    """

    def __init__(self, trajectory_logger: TrajectoryLogger, params: LoggingParams):
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
    def write_interval(self):
        """
        The interval at which log writing is occurring.

        :return: The interval at which log writing is occurring, in steps.
        """
        return self._params.write_interval

    def close(self):
        """
        Close the log.
        """
        self._logger.close()


class ASEOpenMMRunner(NanoverRunner):
    """
    A wrapper class for running an interactive OpenMM simulation with ASE.

    :param simulations OpenMM simulations to run interactively.
    :param params IMD parameters to tune the server.
    :param logging_params Parameters for logging the trajectory of the simulation.
    """

    def __init__(
        self,
        simulations: List[Simulation],
        imd_params: Optional[ImdParams] = None,
        logging_params: Optional[LoggingParams] = None,
    ):
        self._logger = logging.getLogger(__name__)
        self.simulations = simulations
        self.simulation = simulations[0]
        self._validate_simulation()
        if not imd_params:
            imd_params = ImdParams()
        if not logging_params:
            logging_params = LoggingParams()
        self._address = imd_params.address
        self._frame_interval = imd_params.frame_interval
        self._time_step = imd_params.time_step
        self._verbose = imd_params.verbose

        self._sim_index = 0

        self._initialise_calculator(self.simulation, walls=imd_params.walls)
        self._initialise_dynamics()
        self._initialise_server(
            imd_params.name,
            imd_params.address,
            imd_params.port,
            imd_params.discovery,
            imd_params.discovery_port,
        )
        self._initialise_imd(self.app_server, self.dynamics)

        self._initialise_trajectory_logging(logging_params)

        def load(index):
            self._sim_index = int(index % len(self.simulations))
            self.simulation = self.simulations[self._sim_index]
            self._initialise_calculator(self.simulation, walls=imd_params.walls)
            self._initialise_dynamics()
            self.imd.replace_dynamics(self.dynamics)

        def list_():
            return {"simulations": [str(sim) for sim in self.simulations]}

        self.app_server.server.register_command("playback/load", load)
        self.app_server.server.register_command("playback/list", list_)

    @property
    def app_server(self):
        return self._app_server

    @property
    def is_running(self):
        return self.imd.is_running

    def pause(self):
        self.imd.pause()

    def play(self):
        self.imd.play()

    def reset(self):
        self.imd.reset()

    def step(self):
        self.imd.step()

    def _validate_simulation(self):
        """
        Check this runner's simulation for unsupported features and issue the
        relevant warnings.
        """
        if self.simulation.system.getNumConstraints() > 0:
            self._logger.warning(CONSTRAINTS_UNSUPPORTED_MESSAGE)

    @classmethod
    def from_xml(
        cls,
        simulation_xmls: List[str],
        params: Optional[ImdParams] = None,
        logging_params: Optional[LoggingParams] = None,
    ):
        """
        Initialises a :class:`AseOpenMMIMDRunner` from a simulation XML file
        serialised with :func:`serializer.serialize_simulation`.

        :param simulation_xmls: Paths to XML files.
        :param params: The :class: ImdParams to run the server with.
        :param logging_params: The :class:LoggingParams to set up trajectory logging with.
        :return: An OpenMM simulation runner.
        """
        if params is not None:
            platform = params.platform
        else:
            platform = None

        simulations = []
        for path in simulation_xmls:
            with open(path) as infile:
                simulation = serializer.deserialize_simulation(
                    infile.read(),
                    platform_name=platform,
                )
                simulations.append(simulation)

        return cls(simulations, params, logging_params)

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
    def running_discovery(self):
        return self.app_server.running_discovery

    @property
    def discovery_port(self):
        try:
            return self.app_server.discovery.port
        except AttributeError:
            raise AttributeError("Discovery service not running")

    @property
    def dynamics(self):
        """
        Gets the ASE :class:`MolecularDynamics` object that is running the dynamics.

        :return: The ASE molecular dynamics object.
        """
        return self._dynamics

    @property
    def dynamics_interval(self):
        """
        Minimum interval, in seconds,  between frames sent to the frame publisher.
        """
        return self.imd.dynamics_interval

    @dynamics_interval.setter
    def dynamics_interval(self, interval):
        self.imd.dynamics_interval = interval

    def run(
        self,
        steps: Optional[int] = None,
        block: Optional[bool] = None,
        reset_energy: Optional[float] = None,
    ):
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
        self.imd.close()
        if self.logging_info:
            self.logging_info.close()
        self.app_server.close()

    def _initialise_server(
        self,
        name=None,
        address=None,
        port=None,
        run_discovery=True,
        discovery_port=None,
    ):
        address = address or DEFAULT_SERVE_ADDRESS
        if port is None:
            port = DEFAULT_NANOVER_PORT
        server = NanoverServer(address=address, port=port)
        if run_discovery:
            discovery = DiscoveryServer(broadcast_port=discovery_port)
        else:
            discovery = None
        self._app_server = NanoverImdApplication(server, discovery, name)

    def _initialise_imd(self, server, dynamics):
        # set the server to use the OpenMM frame convert for performance purposes.
        self.imd = NanoverASEDynamics(
            server,
            dynamics,
            frame_method=openmm_ase_frame_adaptor,
            frame_interval=self.frame_interval,
        )

    def _initialise_calculator(self, simulation, walls=False):
        self._openmm_calculator = OpenMMCalculator(simulation)
        self._md_calculator = self._openmm_calculator
        self.atoms = self._openmm_calculator.generate_atoms()
        if walls:
            self.atoms.constraints.append(VelocityWallConstraint())
        self.atoms.calc = self._md_calculator

    def _initialise_dynamics(self):
        # Set the momenta corresponding to T=300K
        MaxwellBoltzmannDistribution(self.atoms, temperature_K=300)

        # We do not remove the center of mass (fixcm=False). If the center of
        # mass translations should be removed, then the removal should be added
        # to the OpenMM system.
        self._dynamics = Langevin(
            atoms=self.atoms,
            timestep=self.time_step * units.fs,
            temperature_K=300,
            friction=1e-2,
            fixcm=False,
        )

        if self.verbose:
            self._dynamics.attach(
                MDLogger(
                    self._dynamics,
                    self.atoms,
                    "-",
                    header=True,
                    stress=False,
                    peratom=False,
                ),
                interval=100,
            )

    def _initialise_trajectory_logging(self, logging_params: LoggingParams):
        if not logging_params.trajectory_file:
            self.logging_info = None
            return
        logger = TrajectoryLogger(self.atoms, logging_params.trajectory_file)
        self.imd.on_reset_listeners.append(logger.reset)
        self.logging_info = TrajectoryLoggerInfo(logger, logging_params)
        self.dynamics.attach(logger, logging_params.write_interval)
