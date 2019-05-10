# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.

"""
Interactive molecular dynamics server for use with an ASE molecular dynamics simulation.
"""
from concurrent import futures
from typing import Optional

from ase import Atoms
from ase.calculators.calculator import Calculator
from ase.md.md import MolecularDynamics

from .frame_server import ASEFrameServer
from .imd_calculator import ImdCalculator
from narupa.imd.imd_server import ImdServer
from narupa.trajectory import FrameServer


class ASEImdServer:
    """
    Interactive molecular dynamics runner for use with ASE.

    :param dynamics: A prepared ASE molecular dynamics object to run, with IMD attached.
    :param frame_interval: Interval, in steps, at which to publish frames.
    """

    def __init__(self, dynamics: MolecularDynamics, frame_interval=1):
        self.frame_server = FrameServer()
        self.imd_server = ImdServer()
        self.dynamics = dynamics
        calculator = self.dynamics.atoms.get_calculator()
        self.imd_calculator = ImdCalculator(self.imd_server.service, calculator)
        self.atoms.set_calculator(self.imd_calculator)
        self.dynamics.attach(ASEFrameServer(self.atoms, self.frame_server), interval=frame_interval)
        self.threads = futures.ThreadPoolExecutor(max_workers=1)
        self._run_task = None
        self._cancelled = False

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
        :return:
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
        self._cancelled = True
        if wait:
            self._run_task.result()

    def close(self):
        self.imd_server.close()
        self.frame_server.close()
