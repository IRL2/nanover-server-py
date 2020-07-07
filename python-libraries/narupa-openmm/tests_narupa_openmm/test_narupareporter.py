# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.
"""
Tests for the :class:`narupa.openmm.NarupaReporter`.
"""

# The use of fixture is not recognized by pylint and leads to erroneous warnings.
# pylint: disable=redefined-outer-name

import pytest

from narupa.trajectory import FrameData
from narupa.openmm import NarupaReporter

from .simulation_utils import assert_basic_simulation_topology

# Pylint does not detect the use of the fixture.
from .simulation_utils import basic_simulation  # pylint: disable=unused-import


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
        self.all_sent_frames = []

    def setup_services(self):
        pass

    def send_frame(self, frame_index: int, frame_data: FrameData):
        self.all_sent_frames.append((frame_index, frame_data))


@pytest.mark.parametrize('current_step, report_interval, expected_steps', (
        (0, 1, 1),
        (0, 3, 3),
        (2, 3, 1),
        (3, 3, 3),
))
# describeNextReport is a name that is part of the OpenMM API. It does not
# conform the naming conventions, but it cannot be changed.
def test_describeNextReport(  # pylint: disable=invalid-name
        basic_simulation, current_step, report_interval, expected_steps,
):
    """
    Test that the next report is well described with different report interval
    and different positions in the simulation.
    """
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
    expected_answer = (expected_steps, True, False, False, False, True)
    assert answer == expected_answer


def test_report(basic_simulation):
    """
    Test that two consecutive calls of :meth:`NarupaReporter.report` result in
    the expected calls of the frame server.

    Two calls are required as the first call to the reporter is the only one
    registering the topology. The first and the second calls are then expected
    to behave differently.
    """
    frame_server = MockFrameServer(address='dummy', port=0)
    reporter = NarupaReporter(
        report_interval=1,
        frame_server=frame_server,
    )
    state = basic_simulation.context.getState(getPositions=True)

    reporter.report(basic_simulation, state)
    # The first report sends 2 messages: one with the topology and one with
    # the positions. Subsequent reports only send the positions.
    assert len(frame_server.all_sent_frames) == 2

    frame_index, topology = frame_server.all_sent_frames[0]
    assert frame_index == 0
    assert_basic_simulation_topology(topology)

    reporter.report(basic_simulation, state)
    assert len(frame_server.all_sent_frames) == 3
