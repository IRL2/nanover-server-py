from simtk.openmm.app.topology import Topology
from narupa.protocol.topology.topology_pb2 import TopologyData

from .topology import openmm_topology_to_topology_data
from .frame import openmm_positions_to_frame_data

class NarupaReporter(object):
    """NarupaReporter outputs a series of frames from a Simulation to a narupa server.
    To use it, create a NarupaReporter, then add it to the Simulation's list of reporters.
    """

    _topology : Topology
    _topologyData : TopologyData

    def __init__(self, reportInterval, frameServer):
        self._reportInterval = reportInterval
        self._frameServer = frameServer
        self._topology = None
        self._topologyData = None
        self._frameIndex = 0

    def describeNextReport(self, simulation):
        steps = self._reportInterval - simulation.currentStep%self._reportInterval
        return (steps, True, False, False, False, False)

    def report(self, simulation, state):
        self._topology = simulation.topology
        if(self._frameIndex == 0):
            self._topologyData = openmm_topology_to_topology_data(self._topology)
            self._frameServer.send_topology(self._frameIndex, self._topologyData)
        self._frameData = openmm_positions_to_frame_data(state.getPositions())
        self._frameServer.send_frame(self._frameIndex, self._frameData)
        self._frameIndex += 1