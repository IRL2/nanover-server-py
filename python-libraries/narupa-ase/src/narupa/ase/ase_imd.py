# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.

"""
Interactive molecular dynamics server for use with an ASE molecular dynamics simulation.
"""

from ase import Atoms
from ase.calculators.calculator import Calculator
from ase.md.md import MolecularDynamics

from narupa.ase import NarupaASE
from narupa.ase.imd_calculator import ImdCalculator
from narupa.imd.imd_server import ImdServer
from narupa.trajectory import FrameServer


class ASEImd:
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
        self.dynamics.attach(NarupaASE(self.atoms, self.frame_server), interval=frame_interval)

    @property
    def calculator(self) -> Calculator:
        """
        The calculator being used to compute internal energy and forces.
        :return: ASE calculator.
        """
        return self.imd_calculator.calculator

    @property
    def atoms(self) -> Atoms:
        """
        The atoms in the MD system.
        :return: ASE atoms.
        """
        return self.dynamics.atoms

    @atoms.setter
    def atoms(self, value: Atoms):
        self.dynamics.atoms = value

    def run(self, steps=None):
        """
        Runs the molecular dynamics forward the given number of steps.
        :param steps:
        :return:
        """
        if steps is None:
            while True:
                self.dynamics.run(1000)
        else:
            self.dynamics.run(steps)

    def close(self):
        self.imd_server.close()
        self.frame_server.close()
