from simtk.openmm.app.topology import Topology

from narupa.protocol.topology.topology_pb2 import TopologyData
from narupy.openmm import openmm_to_frame_data


class NarupaReporter(object):
    """NarupaReporter outputs a series of frames from a Simulation to a narupa server.
    To use it, create a NarupaReporter, then add it to the Simulation's list of reporters.
    """

    _topology: Topology
    _topologyData: TopologyData

    def __init__(self, *, report_interval, frame_server):
        self._reportInterval = report_interval
        self._frameServer = frame_server
        self._topology = None
        self._frameData = None
        self._frameIndex = 0

    def describeNextReport(self, simulation):
        steps = self._reportInterval - simulation.currentStep % self._reportInterval
        return steps, True, False, False, False, False

    def report(self, simulation, state):
        self._topology = simulation.topology
        self._frameData = openmm_to_frame_data(positions=state.getPositions(),
                                               topology=self._topology)
        self._frameServer.send_frame(self._frameIndex, self._frameData)
        self._frameIndex += 1
