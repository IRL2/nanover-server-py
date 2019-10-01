# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.

"""
Interactive molecular dynamics server for use with an ASE molecular dynamics simulation.
"""
import logging
from concurrent import futures
from typing import Optional, Callable

from ase import Atoms
from ase.calculators.calculator import Calculator
from ase.md.md import MolecularDynamics

from .frame_server import send_ase_frame
from .imd_calculator import ImdCalculator
from narupa.imd.imd_server import ImdServer
from narupa.trajectory import FrameServer


class ASEImdServer:
    """
    Interactive molecular dynamics runner for use with ASE.

    :param dynamics: A prepared ASE molecular dynamics object to run, with IMD attached.
    :param frame_interval: Interval, in steps, at which to publish frames.
    :param frame_method(ase_atoms, frame_server): Method to use to generate frames, given the the ASE :class: Atoms
           and a :class: FrameServer.
    """

    def __init__(self, dynamics: MolecularDynamics,
                 frame_method:Optional[Callable]=None,
                 frame_interval=1,
                 address:Optional[str]=None,
                 trajectory_port:Optional[int]=None,
                 imd_port:Optional[int]=None):
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

    def run(self, steps: Optional[int] = None):
        """
        Runs the molecular dynamics forward the given number of steps.

        :param steps: If passed, will run the given number of steps, otherwise will run forever
        on a background thread and immediately return.
        """
        if steps is None:
            self._run_task = self.threads.submit(self._run_forever)
        else:
            self.dynamics.run(steps)

    def _run_forever(self):
        while not self._cancelled:
            self.dynamics.run(10)
        self._cancelled = False

    def cancel_run(self, wait=False):
        """
        Cancel molecular dynamics that is running on a background thread.

        :param wait: Whether to block and wait for the molecular dynamics to halt before returning.
        """
        if self._cancelled:
            return
        self._cancelled = True
        if wait:
            self._run_task.result()

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
