# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.
"""
Facilities to run an OpenMM simulation.
"""

import sys

import numpy as np

from simtk.openmm import app

from narupa.openmm import serializer


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
    def from_xml_input(cls, input_xml):
        """
        Create a runner from a serialized simulation.

        :param input_xml: Path to an XML serialised OpenMM simulation.
        :return: An instance of the class.

        .. seealso::

            The XML serialized simulation can be produced by
            :func:`narupa.openmm.serializer.serialize_simulation`.

        """
        with open(str(input_xml)) as infile:
            simulation = serializer.deserialize_simulation(infile.read())
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

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        pass
