# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.
"""
Provide a reporter for OpenMM simulation to publish frames as a Narupa server.
"""
from typing import Union

from simtk.openmm.app.topology import Topology

from narupa.trajectory.frame_server import FrameServer
from .converter import openmm_to_frame_data


class NarupaReporter:
    """
    Outputs a series of frames from a Simulation to a narupa server.

    To use it, create a NarupaReporter, then add it to the Simulation's list
    of reporters.

    Example
    =======

    .. code-block:: python

        frame_server = FrameServer(address="localhost", port=54321)
        frame_reporter = NarupaReporter(report_interval=5,frame_server=frame_server)
        # Assume some OpenMM simulation already exists
        simulation.reporters.add(frame_reporter)

    :param report_interval: Interval in frames between two reports.
    :param frame_server: Instance of a Narupa frame server.
    """
    _topology: Union[Topology, None]

    def __init__(self, *, report_interval, frame_server):
        self._reportInterval = report_interval
        self._frameServer = frame_server
        self._topology = None
        self._frameData = None
        self._frameIndex = 0

    # The name of the method is part of the OpenMM API. It cannot be made to
    # conform PEP8.
    # noinspection PyPep8Naming
    def describeNextReport(self, simulation):  # pylint: disable=invalid-name
        steps = self._reportInterval - simulation.currentStep % self._reportInterval
        # The reporter needs:
        # - the positions
        # - not the velocities
        # - not the forces
        # - not the energies
        # - positions are unwrapped
        return steps, True, False, False, False, True

    def report(self, simulation, state):
        if self._frameIndex == 0:
            self._topology = simulation.topology
            self._frameData = openmm_to_frame_data(state=None,
                                                   topology=self._topology)
            self._frameServer.send_frame(self._frameIndex, self._frameData)
        self._frameData = openmm_to_frame_data(state=state,
                                               topology=None)
        self._frameServer.send_frame(self._frameIndex, self._frameData)
        self._frameIndex += 1
