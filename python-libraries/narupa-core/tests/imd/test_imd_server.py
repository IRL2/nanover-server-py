import time
from unittest.mock import Mock

import grpc
import pytest
from narupa.imd.imd_client import ImdClient, delayed_generator
from narupa.imd.imd_server import ImdServer
from narupa.imd.particle_interaction import ParticleInteraction
from narupa.protocol.imd import InteractionEndReply
from narupa.protocol.imd.imd_pb2_grpc import InteractiveMolecularDynamicsStub


@pytest.fixture
def imd_server():
    server = ImdServer(address=None, port=None)
    yield server
    server.close()


@pytest.fixture
def imd_client():
    client = ImdClient(address='localhost', port=54322)
    yield client
    client.close()


@pytest.fixture
def stub():
    channel = grpc.insecure_channel("{0}:{1}".format('localhost', 54322))
    stub = InteractiveMolecularDynamicsStub(channel)
    yield stub
    channel.close()


@pytest.fixture
def interaction():
    return ParticleInteraction()


@pytest.fixture
def interactions():
    return [ParticleInteraction()] * 10


def test_server(imd_server):
    assert imd_server is not None


def test_publish_interaction(imd_server, stub, interaction):
    mock = Mock()
    imd_server.service.set_callback(mock.callback)
    reply = stub.PublishInteraction((i.proto for i in [interaction]))
    assert isinstance(reply, InteractionEndReply)
    mock.callback.assert_called_once()


def test_publish_multiple_interactions(imd_server, imd_client):
    mock = Mock()
    imd_server.service.set_callback(mock.callback)
    first_set = [ParticleInteraction()] * 10
    second_set = [ParticleInteraction(interaction_id="2")] * 10
    imd_client.publish_interactions_async(delayed_generator(first_set, delay=0.1))
    result = imd_client.publish_interactions(delayed_generator(second_set, delay=0.15))
    assert result is not None
    assert mock.callback.call_count == len(first_set) + len(second_set)


def test_multiplexing_interactions(imd_server, imd_client):
    """
    The server accepts multiplexing interactions, in which different interactions
    are transmitted over the same stream. While not the typical usage, it is tested.
    """
    mock = Mock()
    imd_server.service.set_callback(mock.callback)
    first_set = [ParticleInteraction()] * 10
    second_set = [ParticleInteraction(interaction_id="2")] * 10
    interleaved = [val for pair in zip(first_set, second_set) for val in pair]
    # TODO use a coroutine awaiting input as the generator to control this without needing sleeps
    imd_client.publish_interactions_async(delayed_generator(interleaved, delay=0.01))
    time.sleep(0.04)
    assert len(imd_server.service.interactions) == 2
    time.sleep(0.4)
    assert mock.callback.call_count == len(first_set) + len(second_set)


def test_clear_interactions(imd_server, imd_client, interactions):
    """
    Tests that after interacting the set of interactions are cleared
    """
    mock = Mock()
    imd_server.service.set_callback(mock.callback)
    imd_client.publish_interactions_async(delayed_generator(interactions, delay=0.01))
    time.sleep(0.04)
    assert len(imd_server.service.interactions) == 1
    time.sleep(0.3)
    assert len(imd_server.service.interactions) == 0


def test_repeat_interactions(imd_server, imd_client, interactions):
    mock = Mock()
    imd_server.service.set_callback(mock.callback)
    imd_client.publish_interactions(delayed_generator(interactions, delay=0.01))
    assert mock.callback.call_count == len(interactions)
    imd_client.publish_interactions(delayed_generator(interactions, delay=0.01))
    assert mock.callback.call_count == 2 * len(interactions)


def test_publish_identical_interactions(imd_server, imd_client, interactions):
    """
    Tests that publishing identical interactions at the same time throws a gRPC exception.
    """
    mock = Mock()
    imd_server.service.set_callback(mock.callback)
    imd_client.publish_interactions_async(delayed_generator(interactions, delay=0.1))
    with pytest.raises(grpc.RpcError):
        imd_client.publish_interactions(delayed_generator(interactions, delay=0.15))
