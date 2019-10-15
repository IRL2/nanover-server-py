# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.

"""
Interactive molecular dynamics server for use with an ASE molecular dynamics simulation.
"""
import logging
from concurrent import futures
from typing import Optional, Callable

import numpy as np

from ase import Atoms
from ase.calculators.calculator import Calculator
from ase.lattice.cubic import FaceCenteredCubic
from ase.md import Langevin
from ase.md.md import MolecularDynamics
from narupa.app import NarupaClient

from narupa.ase.frame_server import send_ase_frame
from narupa.ase.imd_calculator import ImdCalculator
from narupa.ase.converter import EV_TO_KJMOL
from narupa.imd.imd_server import ImdServer
from narupa.trajectory import FrameServer


class ASEImdServer:
    """
    Interactive molecular dynamics runner for use with ASE.

    :param dynamics: A prepared ASE molecular dynamics object to run, with IMD attached.
    :param frame_interval: Interval, in steps, at which to publish frames.
    :param frame_method(ase_atoms, frame_server): Method to use to generate frames, given the the ASE :class:`Atoms`
           and a :class: FrameServer.


    Example
    =======

    >>> atoms = FaceCenteredCubic(directions=[[1, 0, 0], [0, 1, 0], [0, 0, 1]], symbol="Cu", size=(2, 2, 2), pbc=True)
    >>> dynamics = Langevin(atoms, timestep=0.5, temperature=300, friction=1.0)
    >>> server = ASEImdServer(dynamics) # create the server with the molecular dynamics object.
    >>> client = NarupaClient(run_multiplayer=False) # have a client connect to the server
    >>> server.run(5) # run some dynamics.
    >>> client.first_frame.particle_count # client will have received some frames!
    32
    >>> # Alternatively, use a 'with' statement to manage the context.
    >>> client.close()
    >>> server.close()

    """

    def __init__(self, dynamics: MolecularDynamics,
                 frame_method: Optional[Callable] = None,
                 frame_interval=1,
                 address: Optional[str] = None,
                 trajectory_port: Optional[int] = None,
                 imd_port: Optional[int] = None):
        if frame_method is None:
            frame_method = send_ase_frame
        self.frame_server = FrameServer(address=address, port=trajectory_port)
        self.imd_server = ImdServer(address=address, port=imd_port)
        self.dynamics = dynamics
        calculator = self.dynamics.atoms.get_calculator()
        self.imd_calculator = ImdCalculator(self.imd_server.service, calculator, dynamics=dynamics)
        self.atoms.set_calculator(self.imd_calculator)
        self.dynamics.attach(frame_method(self.atoms, self.frame_server), interval=frame_interval)
        self.threads = futures.ThreadPoolExecutor(max_workers=1)
        self._run_task = None
        self._cancelled = False

        self._initial_positions = self.atoms.get_positions()
        self._initial_velocities = self.atoms.get_velocities()
        self._initial_box = self.atoms.get_cell()
        self.on_reset_listeners = []

        self.logger = logging.getLogger(__name__)
        self.logger.info(f"Running frame server at {address}:{trajectory_port}")
        self.logger.info(f"Running IMD server at {address}:{imd_port}")

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

    def run(self, steps: Optional[int] = None,
            block: Optional[bool] = None, reset_energy: Optional[float] = None):
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
        # The default is to be blocking if a number of steps is provided, and
        # not blocking if we run forever.
        if block is None:
            block = (steps is not None)
        if steps is None:
            steps = float('inf')
        if block:
            self._run(steps, reset_energy)
        else:
            self._run_task = self.threads.submit(self._run, steps, reset_energy)

    def _run(self, steps, reset_energy):
        remaining_steps = steps
        while not self._cancelled and remaining_steps > 0:
            steps_for_this_iteration = min(10, remaining_steps)
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
        if self._cancelled:
            return
        self._cancelled = True
        if wait:
            self._run_task.result()

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

        """
        self.atoms.set_positions(self._initial_positions)
        self.atoms.set_velocities(self._initial_velocities)
        self.atoms.set_cell(self._initial_box)
        self._call_on_reset()

    def _call_on_reset(self):
        for callback in self.on_reset_listeners:
            callback()

    def close(self):
        """
        Cancels the molecular dynamics if it is running, and closes
        the IMD and frame servers.
        """
        self.cancel_run()
        self.imd_server.close()
        self.frame_server.close()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
