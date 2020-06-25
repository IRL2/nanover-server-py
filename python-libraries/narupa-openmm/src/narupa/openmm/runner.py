# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.
"""
Facilities to run an OpenMM simulation.
"""
from typing import Union, TypeVar, Type, Optional
import sys
import os
import logging

import numpy as np

from simtk.openmm import app

from narupa.openmm import serializer
from narupa.app import NarupaImdApplication
from .imd import NarupaImdReporter, get_imd_forces_from_system, create_imd_force

RunnerClass = TypeVar('RunnerClass', bound='Runner')


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
    def __init__(
            self,
            simulation: app.Simulation,
            name: Optional[str] = None,
            address: Optional[str] = None,
            port: Optional[int] = None,
    ):
        self.simulation = simulation
        self._verbose_reporter = app.StateDataReporter(
            sys.stdout, 10,
            step=True,
            speed=True,
            remainingTime=False,
            potentialEnergy=True,
        )
        self.app = NarupaImdApplication.basic_server(name, address, port)
        putative_imd_forces = get_imd_forces_from_system(simulation.system)
        if not putative_imd_forces:
            raise ValueError(
                'The simulation must include an appropriate force for imd.')
        if len(putative_imd_forces) > 1:
            logging.warning(
                f'More than one force could be used as imd force '
                f'({len(putative_imd_forces)}); taking the last one.'
            )
        # In case there is more than one compatible fore we take the last one.
        # The forces are in the order they have been added, so we take the last
        # one that have been added. This is the most likely to have been added
        # for the purpose of this runner, the other ones are likely leftovers
        # or forces created for another purpose.
        imd_force = putative_imd_forces[-1]
        self.reporter = NarupaImdReporter(
            frame_interval=5,
            force_interval=10,
            imd_force=imd_force,
            imd_service=self.app.imd,
            frame_publisher=self.app.frame_publisher,
        )
        self.simulation.reporters.append(self.reporter)

    @classmethod
    def from_xml_input(
            cls: Type[RunnerClass],
            input_xml: Union[str, bytes, os.PathLike],
            name: Optional[str] = None,
            address: Optional[str] = None,
            port: Optional[int] = None,
    ) -> RunnerClass:
        """
        Create a runner from a serialized simulation.

        :param input_xml: Path to an XML serialised OpenMM simulation.
        :return: An instance of the class.

        .. seealso::

            The XML serialized simulation can be produced by
            :func:`narupa.openmm.serializer.serialize_simulation`.

        """
        imd_force = create_imd_force()
        with open(str(input_xml)) as infile:
            simulation = serializer.deserialize_simulation(
                infile.read(), imd_force=imd_force)
        return cls(simulation)

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

    def run(self, n_steps: float = np.inf) -> None:
        """
        Runs the simulation for a given number of steps, by default an infinity.

        :param n_steps: The number of steps to run. Infinity by default.
        """
        self.simulation.step(n_steps)

    def close(self):
        self.app.close()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
