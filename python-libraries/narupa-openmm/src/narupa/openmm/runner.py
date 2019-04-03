"""
Facilities to run an OpenMM simulation.
"""

import sys

import numpy as np

import simtk.openmm as mm
from simtk.openmm import app

# The prefixed units are programmatically added to the simtk.unit module, thus
# there are not found by pyint or PyCharm.
from simtk.unit import kelvin, picosecond  # pylint: disable=no-name-in-module

# TODO: The Runner class can be constructed from a Simulation instance
#  (__init__) or from a serialized system and a PDB file (from_xml_input).
#  The serialized system does not include the integrator so we need a way to
#  construct the runner from a serialized Simulation.


class Runner:
    """
    Convenience class to run an OpenMM simulation.

    A :class:`Runner` object wraps an OpenMM simulation. The
    :class:`app.Simulation` instance is accessible via the :attr:`simulation`
    attribute.

    Actually starting the simulation is done with the :meth:`run` method. The
    method takes an number of steps to run as an argument; by default, the
    simulation runs indefinitely.

    The verbosity can be adjusted by setting the :attr:`verbose` attribute, or
    by using the :meth:`make_verbose` and :meth:`make_quiet` methods.
    """
    def __init__(self, simulation: app.Simulation):
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
        """
        Construct an instance of the class from a serialized OpenMM system.

        .. warning::

            The serialised OpenMM system does not describe the integrator. As
            a temporary measure, this constructor sets a Langevin integrator at
            300 Kelvin, with an integration step of 2 femtoseconds.

        :param input_xml: Path to an XML serialised OpenMM system.
        :param pdb_path: Path to the corresponding PDB file.
        :return: An instance of the class.
        """
        with open(str(input_xml)) as infile:
            system = mm.XmlSerializer.deserialize(infile.read())
        pdb = app.PDBFile(str(pdb_path))
        integrator = mm.LangevinIntegrator(
            300 * kelvin,
            1 / picosecond,
            0.002 * picosecond
        )
        simulation = app.Simulation(pdb.topology, system, integrator)
        simulation.context.setPositions(pdb.positions)
        return cls(simulation)

    def make_verbose(self):
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

    def make_quiet(self):
        """
        Detach the verbosity reporter if it is attached.

        .. seealso:: :meth:`make_verbose`
        """
        if self.verbose:
            self.simulation.reporters.remove(self._verbose_reporter)

    @property
    def verbose(self):
        """
        Returns ``True`` if the verbosity reporter is attached.
        """
        return self._verbose_reporter in self.simulation.reporters

    @verbose.setter
    def verbose(self, value):
        """
        Sets the verbosity; attach or detach de verbosity reporter if needed.
        """
        if value:
            self.make_verbose()
        else:
            self.make_quiet()

    def run(self, n_steps=np.inf):
        """
        Runs the simulation for a given number of steps, by default an infinity.

        :param n_steps: The number of steps to run. Infinity by default.
        """
        self.simulation.step(n_steps)
