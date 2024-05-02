"""
Facilities to run an OpenMM simulation.
"""

from typing import Union, TypeVar, Type, Optional, Dict
import sys
import os
import logging
import warnings
from concurrent import futures
from threading import RLock
from io import StringIO

import openmm
from openmm import app


from nanover.openmm import serializer
from nanover.app import NanoverImdApplication, NanoverRunner
from .imd import NanoverImdReporter, get_imd_forces_from_system, create_imd_force
from nanover.utilities.event import Event
from nanover.trajectory.frame_server import (
    PLAY_COMMAND_KEY,
    RESET_COMMAND_KEY,
    STEP_COMMAND_KEY,
    PAUSE_COMMAND_KEY,
    GET_DYNAMICS_INTERVAL_COMMAND_KEY,
    SET_DYNAMICS_INTERVAL_COMMAND_KEY,
)
from nanover.utilities.timing import VariableIntervalGenerator

GET_FRAME_INTERVAL_COMMAND_KEY = "trajectory/get-frame-interval"
SET_FRAME_INTERVAL_COMMAND_KEY = "trajectory/set-frame-interval"
GET_FORCE_INTERVAL_COMMAND_KEY = "imd/get-force-interval"
SET_FORCE_INTERVAL_COMMAND_KEY = "imd/set-force-interval"

RunnerClass = TypeVar("RunnerClass", bound="OpenMMRunner")


class OpenMMRunner(NanoverRunner):
    """
    Convenience class to run an OpenMM simulation.

    A :class:`Runner` object wraps an OpenMM simulation. The
    :class:`app.Simulation` instance is accessible via the :attr:`simulation`
    attribute.

    Actually starting the simulation is done with the :meth:`run` method. The
    method takes an number of steps to run as an argument; by default, the
    simulation runs indefinitely.

    By default, one frame is sent to clients every 5 MD steps. This rate can
    be changed by setting the :attr:`frame_interval` attribute. By default, as
    well, the dynamics is throttled to send 30 frames per seconds (therefore
    running 150 MD steps per second). This rate can be changed by setting the
    target interval between frames in seconds with the
    :attr:`dynamics_interval` attribute. Setting the interval to 0 causes the
    dynamics to not be throttled.

    The verbosity can be adjusted by setting the :attr:`verbose` attribute, or
    by using the :meth:`make_verbose` and :meth:`make_quiet` methods.

    :param simulation: The OpenMM simulation to run. It must have an OpenMM
        force object compatible with iMD. This force can be added using
        :func:`nanover.openmm.imd.add_imd_force_to_system` or provided to
        :func:`nanover.openmm.serializer.deserialize_simulation` with the
        ``imd_force`` argument.
    :param name: A friendly name for the runner. It will be displayed by ESSD.
    :param address: The IP address the server binds to.
    :param port: The port the server listens to.
    """

    def __init__(
        self,
        simulation: app.Simulation,
        name: Optional[str] = None,
        address: Optional[str] = None,
        port: Optional[int] = None,
    ):
        self.simulation = simulation
        self._verbose_reporter = app.StateDataReporter(
            sys.stdout,
            10,
            step=True,
            speed=True,
            remainingTime=False,
            potentialEnergy=True,
        )
        self._app_server = NanoverImdApplication.basic_server(name, address, port)
        potential_imd_forces = get_imd_forces_from_system(simulation.system)
        if not potential_imd_forces:
            raise ValueError(
                "The simulation must include an appropriate force for imd."
            )
        if len(potential_imd_forces) > 1:
            logging.warning(
                f"More than one force could be used as imd force "
                f"({len(potential_imd_forces)}); taking the last one."
            )
        # In case there is more than one compatible force we take the last one.
        # The forces are in the order they have been added, so we take the last
        # one that have been added. This is the most likely to have been added
        # for the purpose of this runner, the other ones are likely leftovers
        # or forces created for another purpose.
        imd_force = potential_imd_forces[-1]
        self.reporter = NanoverImdReporter(
            frame_interval=5,
            force_interval=10,
            imd_force=imd_force,
            imd_state=self.app_server.imd,
            frame_publisher=self.app_server.frame_publisher,
        )
        self.simulation.reporters.append(self.reporter)

        initial_state_fake_file = StringIO()
        self.simulation.saveState(initial_state_fake_file)
        self._initial_state = initial_state_fake_file.getvalue()

        self._variable_interval_generator = VariableIntervalGenerator(1 / 30)
        self.threads = futures.ThreadPoolExecutor(max_workers=1)
        self._cancel_lock = RLock()
        self._run_task: Optional[futures.Future[None]] = None
        self._cancelled = False

        self.on_reset = Event()

        self._register_commands()

    @classmethod
    def from_xml_input(
        cls: Type[RunnerClass],
        input_xml: Union[str, bytes, os.PathLike],
        name: Optional[str] = None,
        address: Optional[str] = None,
        port: Optional[int] = None,
        platform: Optional[str] = None,
    ) -> RunnerClass:
        """
        Create a runner from a serialized simulation.

        :param input_xml: Path to an XML serialised OpenMM simulation.
        :param name: A friendly name for the runner. It will be displayed
            by ESSD.
        :param address: The IP address the server binds to.
        :param port: The port the server listens to.
        :return: An instance of the class.

        .. seealso::

            The XML serialized simulation can be produced by
            :func:`nanover.openmm.serializer.serialize_simulation`.

        """
        imd_force = create_imd_force()
        with open(str(input_xml)) as infile:
            simulation = serializer.deserialize_simulation(
                infile.read(), imd_force=imd_force, platform_name=platform
            )
        return cls(simulation, name=name, address=address, port=port)

    @property
    def app_server(self):
        return self._app_server

    @property
    def app(self):
        warnings.warn(
            'The property "app" is deprecated and will be removed in '
            'a later version. Use "app_server" instead.',
            DeprecationWarning,
        )
        return self.app_server

    @property
    def frame_interval(self) -> int:
        """
        Send a frame every N steps.
        """
        return self.reporter.frame_interval

    @frame_interval.setter
    def frame_interval(self, interval: int) -> None:
        self.reporter.frame_interval = interval

    @property
    def force_interval(self) -> int:
        """
        Update iMD interactions every N steps.
        """
        return self.reporter.force_interval

    @force_interval.setter
    def force_interval(self, interval: int) -> None:
        self.reporter.force_interval = interval

    @property
    def verbosity_interval(self) -> int:
        """
        Display the verbosity report every N steps.

        If the runner is not verbose, then the verbosity interval is 0.
        Same wise, if this interval is set to 0, then the runner is made quiet.
        """
        if self.verbose:
            return self._verbose_reporter._reportInterval
        return 0

    @verbosity_interval.setter
    def verbosity_interval(self, interval: int) -> None:
        if interval:
            self._verbose_reporter._reportInterval = interval
            self.make_verbose()
        else:
            self.make_quiet()

    def make_verbose(self) -> None:
        """
        Attach a verbosity reporter if it is not already attached.

        The verbosity reporter reports the step number, the potential energy
        in kJ/mol, and the simulation speed in ns/day. This report is displayed
        every 10 simulation steps.

        .. seealso::

            The :meth:`make_quiet` method removes the verbosity reporter.

        """
        if not self.verbose:
            self.simulation.reporters.append(self._verbose_reporter)

    def make_quiet(self) -> None:
        """
        Detach the verbosity reporter if it is attached.

        .. seealso:: :meth:`make_verbose`
        """
        if self.verbose:
            self.simulation.reporters.remove(self._verbose_reporter)

    @property
    def verbose(self) -> bool:
        """
        Returns ``True`` if the verbosity reporter is attached.
        """
        return self._verbose_reporter in self.simulation.reporters

    @verbose.setter
    def verbose(self, value: bool):
        """
        Sets the verbosity; attach or detach the verbosity reporter if needed.
        """
        if value:
            self.make_verbose()
        else:
            self.make_quiet()

    @property
    def is_running(self) -> bool:
        """
        Whether or not the molecular dynamics is currently running on a
        background thread or not.
        :return: `True`, if molecular dynamics is running, `False` otherwise.
        """
        # ideally we'd just check _run_task.running(), but there can be a delay
        # between the task starting and hitting the running state.
        return self._run_task is not None and not (
            self._run_task.cancelled() or self._run_task.done()
        )

    @property
    def dynamics_interval(self):
        """
        Minimum interval, in seconds,  between frames sent to the frame publisher.
        """
        return self._variable_interval_generator.interval

    @dynamics_interval.setter
    def dynamics_interval(self, interval):
        self._variable_interval_generator.interval = interval

    def run(
        self,
        steps: Optional[int] = None,
        block: Optional[bool] = None,
    ) -> None:
        if self.is_running:
            raise RuntimeError("Dynamics are already running on a thread!")
        # The default is to be blocking if a number of steps is provided, and
        # not blocking if we run forever.
        if block is None:
            block = steps is not None
        if block:
            self._run(steps)
        else:
            self._run_task = self.threads.submit(self._run, steps)

    def _run(self, steps: Optional[int]) -> None:
        try:
            remaining_steps = steps if steps is not None else float("inf")
            for _ in self._variable_interval_generator.yield_interval():
                if self._cancelled or remaining_steps <= 0:
                    break
                steps_for_this_iteration = min(self.frame_interval, remaining_steps)
                try:
                    self.simulation.step(steps_for_this_iteration)
                except (ValueError, openmm.OpenMMException):
                    # We want to stop running if the simulation exploded in a way
                    # that prevents OpenMM to run. Otherwise, we will be a a state
                    # where OpenMM raises an exception which would make the runner
                    # unusable. The OpenMMException is typically raised by OpenMM
                    # itself when something is NaN; the ValueError is typically
                    # raised by the StateReporter when the energy is NaN.
                    break
                remaining_steps -= steps_for_this_iteration
            self._cancelled = False
        except Exception as err:
            print(f"Error whith run: {err}")

    def step(self):
        with self._cancel_lock:
            self.cancel_run(wait=True)
            self.run(self.frame_interval, block=True)
            self.cancel_run(wait=True)

    def pause(self):
        with self._cancel_lock:
            self.cancel_run(wait=True)

    def play(self):
        with self._cancel_lock:
            self.cancel_run(wait=True)
        self.run()

    def reset(self):
        with self._cancel_lock:
            was_running = self.is_running
            self.cancel_run(wait=True)
            self.simulation.context.reinitialize()
            initial_state_fake_file = StringIO(self._initial_state)
            self.simulation.loadState(initial_state_fake_file)
            self.on_reset.invoke()
        if was_running:
            self.run()
        else:
            self.step()

    def cancel_run(self, wait: bool = False) -> None:
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

    def _set_frame_interval(self, interval: int) -> None:
        self.frame_interval = int(interval)

    def _get_frame_interval(self) -> Dict[str, int]:
        return {"interval": self.frame_interval}

    def _set_force_interval(self, interval: int) -> None:
        self.force_interval = int(interval)

    def _get_force_interval(self) -> Dict[str, int]:
        return {"interval": self.force_interval}

    def _set_dynamics_interval(self, interval: float) -> None:
        self.dynamics_interval = float(interval)

    def _get_dynamics_interval(self) -> Dict[str, float]:
        return {"interval": self.dynamics_interval}

    def close(self):
        self.cancel_run()
        self.app_server.close()

    def _register_commands(self):
        server = self.app_server.server
        server.register_command(PLAY_COMMAND_KEY, self.play)
        server.register_command(RESET_COMMAND_KEY, self.reset)
        server.register_command(STEP_COMMAND_KEY, self.step)
        server.register_command(PAUSE_COMMAND_KEY, self.pause)
        server.register_command(
            SET_FRAME_INTERVAL_COMMAND_KEY, self._set_frame_interval
        )
        server.register_command(
            GET_FRAME_INTERVAL_COMMAND_KEY, self._get_frame_interval
        )
        server.register_command(
            SET_FORCE_INTERVAL_COMMAND_KEY, self._set_force_interval
        )
        server.register_command(
            GET_FORCE_INTERVAL_COMMAND_KEY, self._get_force_interval
        )
        server.register_command(
            SET_DYNAMICS_INTERVAL_COMMAND_KEY, self._set_dynamics_interval
        )
        server.register_command(
            GET_DYNAMICS_INTERVAL_COMMAND_KEY, self._get_dynamics_interval
        )
