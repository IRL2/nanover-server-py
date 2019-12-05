import time

import pytest
from typing import NamedTuple, Sequence, Callable
from ..multiplayer import server_client_pair
from narupa.core.grpc_utils import subscribe_channel_connectivity_change
from grpc import insecure_channel, ChannelConnectivity

NOMINAL_WAIT_TIME = 0.01


class ConnectivityRecorder(NamedTuple):
    history: Sequence[ChannelConnectivity]
    callback: Callable[[ChannelConnectivity], None]


@pytest.fixture
def connectivity_recorder() -> ConnectivityRecorder:
    connectivity_history = []

    def record_connectivity(connectivity: ChannelConnectivity):
        print(connectivity)
        connectivity_history.append(connectivity)

    return ConnectivityRecorder(connectivity_history, record_connectivity)


def test_subscribe_unused_channel(connectivity_recorder):
    """
    Test that subscribing a unused channel's connectivity calls the callback
    with the default `IDLE` state.
    """
    channel = insecure_channel("localhost:0")
    subscribe_channel_connectivity_change(channel,
                                          connectivity_recorder.callback)
    time.sleep(NOMINAL_WAIT_TIME)
    assert connectivity_recorder.history == [ChannelConnectivity.IDLE]


def test_subscribe_unused_channel_closed(connectivity_recorder):
    """
    Test that subscribing a closed channel's connectivity does not call the
    callback. This will cause an exception on another thread, but it can't be
    tested for.
    """
    channel = insecure_channel("localhost:0")
    channel.close()
    subscribe_channel_connectivity_change(channel,
                                          connectivity_recorder.callback)
    time.sleep(NOMINAL_WAIT_TIME)
    assert connectivity_recorder.history == []


def test_subscribe_unused_channel_connect(connectivity_recorder):
    """
    Test that subscribing a unused channel's connectivity and forcing connection
    when no corresponding server exists calls the callback with `IDLE`,
     `CONNECTING` then `TRANSIENT_FAILURE` states.
    """
    channel = insecure_channel("localhost:0")
    subscribe_channel_connectivity_change(channel,
                                          connectivity_recorder.callback,
                                          force_connection=True)
    time.sleep(NOMINAL_WAIT_TIME)
    assert connectivity_recorder.history == [
        ChannelConnectivity.IDLE,
        ChannelConnectivity.CONNECTING,
        ChannelConnectivity.TRANSIENT_FAILURE,
    ]
