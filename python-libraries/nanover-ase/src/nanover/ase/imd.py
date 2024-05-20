"""
Interactive molecular dynamics server for use with an ASE molecular dynamics simulation.
"""

import logging
from concurrent import futures
from contextlib import contextmanager
from threading import RLock
from typing import Optional, Callable, List, Dict, Any

import numpy as np
from ase import Atoms  # type: ignore
from ase.calculators.calculator import Calculator
from ase.md.md import MolecularDynamics
from nanover.app import NanoverImdApplication
from nanover.ase.converter import EV_TO_KJMOL
from nanover.ase.frame_adaptor import send_ase_frame
from nanover.ase.imd_calculator import ImdCalculator
from nanover.trajectory.frame_server import (
    PLAY_COMMAND_KEY,
    RESET_COMMAND_KEY,
    STEP_COMMAND_KEY,
    PAUSE_COMMAND_KEY,
    GET_DYNAMICS_INTERVAL_COMMAND_KEY,
    SET_DYNAMICS_INTERVAL_COMMAND_KEY,
)
from nanover.utilities.timing import VariableIntervalGenerator


class NanoverASEDynamics:
    """
    Interactive molecular dynamics adaptor for use with ASE.

    :param nanover_imd_app: A :class:`NanoverImdApplication` to pass frames to and
        read forces from.
    :param dynamics: A prepared ASE molecular dynamics object to run,
        with IMD attached.
    :param frame_interval: Interval, in steps, at which to publish frames.
    :param frame_method: Method to use to generate frames, given the ASE
        :class:`Atoms` and a :class:`FramePublisher`. The signature of the
        callback is expected to be ``frame_method(ase_atoms, frame_publisher)``.


    Example
    =======

    >>> from ase.calculators.emt import EMT
    >>> from ase.lattice.cubic import FaceCenteredCubic
    >>> atoms = FaceCenteredCubic(directions=[[1, 0, 0], [0, 1, 0], [0, 0, 1]], symbol="Cu", size=(2, 2, 2), pbc=True)
    >>> atoms.calc = EMT()
    >>> ase_dynamics = Langevin(atoms, timestep=0.5, temperature_K=300, friction=1.0)
    >>> with NanoverASEDynamics.basic_imd(ase_dynamics) as sim: # run basic NanoVer server
    ...
    ...     with NanoverImdClient.autoconnect() as client: # connect an iMD client.
    ...         sim.run(10) # run some dynamics
    ...         client.first_frame.particle_count # the client will have received some MD data!
    32

    """

    on_reset_listeners: List[Callable[[], None]]
    _run_task: Optional[futures.Future]

    dynamics: MolecularDynamics
    imd_calculator: ImdCalculator
    _initial_positions: Any
    _initial_velocities: Any
    _initial_box: Any

    def __init__(
        self,
        nanover_imd_app: NanoverImdApplication,
        dynamics: MolecularDynamics,
        frame_method: Optional[Callable] = None,
        frame_interval=1,
    ):
        if frame_method is None:
            frame_method = send_ase_frame

        self._frame_method = frame_method
        self._imd = nanover_imd_app.imd
        self._server = nanover_imd_app.server
        self._frame_publisher = nanover_imd_app.frame_publisher

        self._cancel_lock = RLock()
        self._register_commands()

        self._variable_interval_generator = VariableIntervalGenerator(1 / 30)
        self._frame_interval = frame_interval
        self.threads = futures.ThreadPoolExecutor(max_workers=1)
        self._run_task = None
        self._cancelled = False

        self.on_reset_listeners = []
        self.logger = logging.getLogger(__name__)

        self.replace_dynamics(dynamics, reset=False)

    @classmethod
    @contextmanager
    def basic_imd(
        cls,
        dynamics: MolecularDynamics,
        address: Optional[str] = None,
        port: Optional[int] = None,
        **kwargs,
    ):
        """
        Initialises basic interactive molecular dynamics running a NanoVer server
        at the given address and port.

        :param dynamics: Molecular dynamics object to attach the server to.
        :param address: Address to run the server at.
        :param port: Port to run the server on.
        :param kwargs: Key-word arguments to pass to the constructor of :class:~NanoverASEDynamics
        :return: Instantiation of a :class:~NanoverASEDynamics configured with the given server parameters and dynamics.
        """
        with NanoverImdApplication.basic_server(address=address, port=port) as app:
            with cls(app, dynamics, **kwargs) as imd:
                yield imd

    @property
    def address(self) -> str:
        """
        The address of the NanoVer server.
        """
        return self._server.address

    @property
    def port(self) -> str:
        """
        The port of the NanoVer server.
        :return:
        """
        return self._server.port

    @property
    def internal_calculator(self) -> Calculator:
        """
        The internal calculator being used to compute internal energy and forces.

        :return: ASE internal calculator.
        """
        return self.imd_calculator.calculator

    @property
    def atoms(self) -> Atoms:
        """
        The atoms in the MD system.

        :return: ASE atoms.
        """
        return self.dynamics.atoms

    @property
    def dynamics_interval(self):
        """
        Minimum interval, in seconds,  between frames sent to the frame publisher.
        """
        return self._variable_interval_generator.interval

    @dynamics_interval.setter
    def dynamics_interval(self, interval):
        self._variable_interval_generator.interval = interval

    @property
    def is_running(self):
        """
        Whether or not the molecular dynamics is currently running on a background thread or not.
        :return: `True`, if molecular dynamics is running, `False` otherwise.
        """
        # ideally we'd just check _run_task.running(), but there can be a delay between the task
        # starting and hitting the running state.
        return self._run_task is not None and not (
            self._run_task.cancelled() or self._run_task.done()
        )

    def step(self):
        """
        Take a single step of the simulation and stop.

        This method is called whenever a client runs the step command, described in :mod:nanover.trajectory.frame_server.
        """
        with self._cancel_lock:
            self.cancel_run(wait=True)
            self.run(self._frame_interval, block=True)
            self.cancel_run(wait=True)

    def pause(self):
        """
        Pause the simulation, by cancelling any current run.

        This method is called whenever a client runs the pause command,
        described in :mod:nanover.trajectory.frame_server.
        """
        with self._cancel_lock:
            self.cancel_run(wait=True)

    def play(self):
        """
        Run the simulation indefinitely

        Cancels any current run and then begins running the simulation on a background thread.

        This method is called whenever a client runs the play command,
        described in :mod:nanover.trajectory.frame_server.

        """
        with self._cancel_lock:
            self.cancel_run(wait=True)
        self.run()

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
        :param block: If ``False``, run in a separate thread. By default, "block"
            is ``None``, which means it is automatically set to ``True`` if a
            number of steps is provided and to ``False`` otherwise.
        :param reset_energy: Threshold of total energy in kJ/mol above which
            the simulation is reset to its initial conditions. If a value is
            provided, the simulation is reset if the total energy is greater
            than this value, or if the total energy is `nan` or infinite. If
            ``None`` is provided instead, then the simulation will not be
            automatically reset.
        """
        if self.is_running:
            raise RuntimeError("Dynamics are already running on a thread!")
        # The default is to be blocking if a number of steps is provided, and
        # not blocking if we run forever.
        if block is None:
            block = steps is not None
        if block:
            self._run(steps, reset_energy)
        else:
            self._run_task = self.threads.submit(self._run, steps, reset_energy)

    def _run(self, steps: Optional[int], reset_energy: Optional[float]):
        remaining_steps = steps if steps is not None else float("inf")
        for _ in self._variable_interval_generator.yield_interval():
            if self._cancelled or remaining_steps <= 0:
                break
            steps_for_this_iteration = min(self._frame_interval, remaining_steps)
            self.dynamics.run(steps_for_this_iteration)
            remaining_steps -= steps_for_this_iteration
            self._reset_if_required(reset_energy)
        self._cancelled = False

    def _reset_if_required(self, reset_energy):
        if reset_energy is not None:
            energy = self.dynamics.atoms.get_total_energy() * EV_TO_KJMOL
            if not np.isfinite(energy) or energy > reset_energy:
                self.reset()

    def cancel_run(self, wait=False):
        """
        Cancel molecular dynamics that is running on a background thread.

        :param wait: Whether to block and wait for the molecular dynamics to
            halt before returning.
        """
        if self._run_task is None:
            return

        if self._cancelled:
            return
        self._cancelled = True
        if wait:
            self._run_task.result()
            self._cancelled = False

    def reset(self):
        """
        Reset the positions, velocities, and box to their initial values.

        When this happens, the "on_reset" event is triggered and all the
        callbacks listed in the :attr:`on_reset_listeners` are called. These
        callbacks are called without arguments and no return value is stored.

        .. note::

            Only the positions, the velocities, and the simulation box are
            reset to their initial values. If a simulation needs any other
            state to be reset or updated, one should register a callback in
            the :attr:`on_reset_listeners` list. The callbacks are executed
            in the order of the list, after the positions, velocities, and box
            are reset.

            Such callbacks also allow to modify the simulation at each reset.
            They would allow, for instance, to draw new velocities, or to
            place molecules differently.

        This method is called whenever a client runs the reset command,
        described in :mod:`nanover.trajectory.frame_server`.

        """
        self.atoms.set_positions(self._initial_positions)
        self.atoms.set_velocities(self._initial_velocities)
        self.atoms.set_cell(self._initial_box)
        self._call_on_reset()

    def replace_dynamics(self, dynamics: MolecularDynamics, reset=True):
        if reset:
            self.reset()

        with self._cancel_lock:
            self.cancel_run(wait=True)

        # replace dynamics and calculator
        self.dynamics = dynamics
        self.imd_calculator = ImdCalculator(
            self._imd,
            self.dynamics.atoms.calc,
            dynamics=dynamics,
        )
        self.atoms.calc = self.imd_calculator

        # replace previous frame method with fresh instance
        self.dynamics.observers.clear()
        self.dynamics.attach(
            self._frame_method(self.atoms, self._frame_publisher),
            interval=self._frame_interval,
        )

        # backup initial state for resets
        self._initial_positions = self.atoms.get_positions()
        self._initial_velocities = self.atoms.get_velocities()
        self._initial_box = self.atoms.get_cell()

    def _call_on_reset(self):
        for callback in self.on_reset_listeners:
            callback()

    def _set_dynamics_interval(self, interval: float) -> None:
        self.dynamics_interval = float(interval)

    def _get_dynamics_interval(self) -> Dict[str, float]:
        return {"interval": self.dynamics_interval}

    def _register_commands(self):
        self._server.register_command(PLAY_COMMAND_KEY, self.play)
        self._server.register_command(RESET_COMMAND_KEY, self.reset)
        self._server.register_command(STEP_COMMAND_KEY, self.step)
        self._server.register_command(PAUSE_COMMAND_KEY, self.pause)
        # TODO: Fix type annotations for commands with arguments
        self._server.register_command(SET_DYNAMICS_INTERVAL_COMMAND_KEY, self._set_dynamics_interval)  # type: ignore
        self._server.register_command(GET_DYNAMICS_INTERVAL_COMMAND_KEY, self._get_dynamics_interval)  # type: ignore

    def close(self):
        """
        Cancels the molecular dynamics if it is running.
        """
        self.cancel_run()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
