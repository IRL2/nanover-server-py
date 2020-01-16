import time

import pytest
from typing import NamedTuple, Sequence, Callable
from narupa.utilities.grpc_utilities import subscribe_channel_connectivity_change
from grpc import insecure_channel, ChannelConnectivity

NOMINAL_WAIT_TIME = 0.01


class ConnectivityRecorder(NamedTuple):
    history: Sequence[ChannelConnectivity]
    callback: Callable[[ChannelConnectivity], None]


@pytest.fixture
def channel():
    with insecure_channel('localhost:0') as channel:
        yield channel


@pytest.fixture
def connectivity_recorder() -> ConnectivityRecorder:
    connectivity_history = []

    def record_connectivity(connectivity: ChannelConnectivity):
        connectivity_history.append(connectivity)

    return ConnectivityRecorder(connectivity_history, record_connectivity)


def test_subscribe_unused_channel(channel, connectivity_recorder):
    """
    Test that subscribing a unused channel's connectivity calls the callback
    with the default `IDLE` state.
    """
    subscribe_channel_connectivity_change(channel,
                                          connectivity_recorder.callback)
    time.sleep(NOMINAL_WAIT_TIME)
    assert connectivity_recorder.history == [ChannelConnectivity.IDLE]


def test_subscribe_unused_channel_closed(channel, connectivity_recorder):
    """
    Test that subscribing a closed channel's connectivity does not call the
    callback. This will cause an exception on another thread, but it can't be
    tested for.
    """
    channel.close()
    subscribe_channel_connectivity_change(channel,
                                          connectivity_recorder.callback)
    time.sleep(NOMINAL_WAIT_TIME)
    assert connectivity_recorder.history == []


def test_subscribe_unused_channel_connect(channel, connectivity_recorder):
    """
    Test that subscribing a unused channel's connectivity and forcing connection
    when no corresponding server exists calls the callback first with `IDLE` and
    ends with a `TRANSIENT_FAILURE` state, potentially entering `CONNECTING`
    state between them.
    """
    subscribe_channel_connectivity_change(channel,
                                          connectivity_recorder.callback,
                                          force_connection=True)
    time.sleep(NOMINAL_WAIT_TIME)
    assert connectivity_recorder.history[0] == ChannelConnectivity.IDLE
    assert connectivity_recorder.history[-1] == ChannelConnectivity.TRANSIENT_FAILURE
