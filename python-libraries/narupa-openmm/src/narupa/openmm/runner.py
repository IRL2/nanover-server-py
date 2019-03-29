import sys

import numpy as np

import simtk.openmm as mm
from simtk.openmm import app
from simtk.unit import kelvin, picosecond


class Runner:
    def __init__(self, simulation):
        self.simulation = simulation
        self._verbose_reporter = app.StateDataReporter(
            sys.stdout, 10,
            step=True,
            speed=True,
            remainingTime=False,
            potentialEnergy=True,
        )

    @classmethod
    def from_xml_input(cls, input_xml, pdb_path):
        with open(input_xml) as infile:
            system = mm.XmlSerializer.deserialize(infile.read())
        pdb = app.PDBFile(pdb_path)
        integrator = mm.LangevinIntegrator(
            300 * kelvin,
            1 / picosecond,
            0.002 * picosecond
        )
        simulation = app.Simulation(pdb.topology, system, integrator)
        simulation.context.setPositions(pdb.positions)
        return cls(simulation)

    def make_verbose(self):
        if not self.verbose:
            self.simulation.reporters.append(self._verbose_reporter)

    def make_quiet(self):
        if self.verbose:
            self.simulation.reporters.remove(self._verbose_reporter)

    @property
    def verbose(self):
        return self._verbose_reporter in self.simulation.reporters

    @verbose.setter
    def verbose(self, value):
        if value:
            self.make_verbose()
        else:
            self.make_quiet()

    def run(self, n_steps=np.inf):
        self.simulation.step(n_steps)
