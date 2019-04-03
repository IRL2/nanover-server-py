"""
Tests for the :class:`narupa.openmm.NarupaReporter`.
"""

import pytest

from narupa.protocol.trajectory.frame_pb2 import FrameData
from narupa.protocol.topology.topology_pb2 import TopologyData

from narupa.openmm import NarupaReporter

from simulation_utils import basic_simulation


class MockFrameServer:
    """
    Pretend to be a :class:`narupa.trajectory.FrameServer`.

    Mocks a :class:`narupa.trajectory.FrameServer`. Instead of sending data to
    a client, this "server" only stores the frame and topologies supposed that
    should be sent.

    The data supposed to be sent is accessible from the
    :attr:`all_sent_topologies` and the :attr:`all_sent_frames` attributes
    as a list of (frame index, :class:`TopologyData`) tuples and a list of
    (frame index, :class:`FrameData`), respectively. In both cases, the frame
    index and the data are the argument passed by the caller.

    The initialisation expects a host address and a host port, which are stored
    under the :attr:`address` and the :attr:`port` attributes, respectively.
    These arguments mimic the signature of the original frame server; they are
    so far to spot changes in this signature.
    """
    def __init__(self, *, address: str, port: int):
        self.address = address
        self.port = port
        self.all_sent_topologies = []
        self.all_sent_frames = []

    def setup_services(self):
        pass

    def send_topology(self, frame_index: int, topology_data: TopologyData):
        self.all_sent_topologies.append((frame_index, topology_data))

    def send_frame(self, frame_index: int, frame_data: FrameData):
        self.all_sent_frames.append((frame_index, frame_data))


@pytest.mark.parametrize('current_step, report_interval, expected_steps', (
        (0, 1, 1),
        (0, 3, 3),
        (2, 3, 1),
        (3, 3, 3),
))
def test_describeNextReport(  # pytest: disable=invalid-name
        basic_simulation, current_step, report_interval, expected_steps,
):
    frame_server = MockFrameServer(address='dummy', port=0)
    reporter = NarupaReporter(
        report_interval=report_interval,
        frame_server=frame_server,
    )
    basic_simulation.currentStep = current_step
    answer = reporter.describeNextReport(basic_simulation)
    # The booleans in the answer are whether or not the reported needs:
    # * positions,
    # * velocities,
    # * forces,
    # * or energies
    # to be pulled from the context. The NarupaReported only needs the
    # positions.
    expected_answer = (expected_steps, True, False, False, False, False)
    assert answer == expected_answer


@pytest.mark.skip(reason='Test system is not good enough, yet.')
def test_report(basic_simulation):
    frame_server = MockFrameServer(address='dummy', port=0)
    reporter = NarupaReporter(
        report_interval=1,
        frame_server=frame_server,
    )
    state = basic_simulation.context.getState(getPositions=True)

    reporter.report(basic_simulation, state)
    assert len(frame_server.all_sent_topologies) == 1
    assert len(frame_server.all_sent_frames) == 1

    topology = frame_server.all_sent_topologies[0]
    # TODO: we need a more complex test system with bonds and more
    #  than one residue
    assert topology.arrays['residue.id'].string_values.values == ['RES']
    assert topology.arrays['residue.chain'].index_values.values == [0]
    assert topology.arrays['atom.id'].index_values.values == [0, 1]
    assert topology.arrays['atom.element'].index_values.values == [15, 16]
    assert topology.arrays['atom.residue'].index_values.values == [0, 0]
    assert topology.arrays['bond'].index_values.values == []

