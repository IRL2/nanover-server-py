"""
Provide a reporter for OpenMM simulation to publish frames as a Narupa server.
"""
from simtk.openmm.app.topology import Topology

from narupa.trajectory.frame_server import FrameServer
from .converter import openmm_to_frame_data


class NarupaReporter:
    """
    Outputs a series of frames from a Simulation to a narupa server.

    To use it, create a NarupaReporter, then add it to the Simulation's list
    of reporters.

    :param report_interval: Interval in frames between two reports.
    :param frame_server: Instance of a Narupa frame server.
    """
    _topology: Topology

    def __init__(self, *, report_interval, frame_server):
        self._reportInterval = report_interval
        self._frameServer = frame_server
        self._topology = None
        self._frameData = None
        self._frameIndex = 0

    # The name of the method is part of the OpenMM API. It cannot be made to
    # conform PEP8.
    def describeNextReport(self, simulation):  # pylint: disable=invalid-name
        steps = self._reportInterval - simulation.currentStep % self._reportInterval
        return steps, True, False, False, False, False

    def report(self, simulation, state):
        if self._frameIndex == 0:
            self._topology = simulation.topology
            self._frameData = openmm_to_frame_data(positions=None,
                                                   topology=self._topology)
            self._frameServer.send_frame(self._frameIndex, self._frameData)
        self._frameData = openmm_to_frame_data(positions=state.getPositions(),
                                               topology=None)
        self._frameServer.send_frame(self._frameIndex, self._frameData)
        self._frameIndex += 1
