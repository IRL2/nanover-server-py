"""
Provide a reporter for OpenMM simulation to publish frames as a Narupa server.
"""
from simtk.openmm.app.topology import Topology
from narupa.protocol.topology.topology_pb2 import TopologyData
from narupa.trajectory.frame_server import FrameServer

from .converters import (
    openmm_positions_to_frame_data,
    openmm_topology_to_topology_data,
)


class NarupaReporter:
    """
    Outputs a series of frames from a Simulation to a narupa server.

    To use it, create a NarupaReporter, then add it to the Simulation's list
    of reporters.

    :param report_interval: Interval in frames between two reports.
    :param frame_server: Instance of a Narupa frame server.
    """

    _topology: Topology
    _topology_data: TopologyData

    def __init__(self, report_interval: int, frame_server: FrameServer):
        self._report_interval = report_interval
        self._frame_server = frame_server
        self._topology = None
        self._topology_data = None
        self._frame_index = 0
        self._frame_data = None

    # The name of the method is part of the OpenMM API. It cannot be made to
    # conform PEP8.
    def describeNextReport(self, simulation):  # pylint: disable=invalid-name
        steps = self._report_interval - simulation.currentStep % self._report_interval
        return steps, True, False, False, False, False

    def report(self, simulation, state):
        self._topology = simulation.topology
        if self._frame_index == 0:
            self._topology_data = openmm_topology_to_topology_data(self._topology)
            self._frame_server.send_topology(self._frame_index, self._topology_data)
        self._frame_data = openmm_positions_to_frame_data(state.getPositions())
        self._frame_server.send_frame(self._frame_index, self._frame_data)
        self._frame_index += 1
