import time

import grpc
import pytest
from .test_frame_server import simple_frame_data, frame_server
from .imd.test_imd_server import imd_client, imd_server, interaction
from narupa.app.client import NarupaClient
import numpy as np


@pytest.fixture
def client():
    client = NarupaClient()
    yield client
    client.close()


@pytest.fixture
def client_server(frame_server, imd_server):
    client = NarupaClient()
    yield client, frame_server, imd_server
    client.close()


def test_receive_frames(client_server, simple_frame_data):
    client, frame_server, imd_server = client_server
    time.sleep(0.2)
    frame_server.send_frame(0, simple_frame_data)
    time.sleep(0.5)
    assert client.first_frame is not None
    assert client.latest_frame == client.first_frame
    assert len(client.frames) == 1


def test_receive_multiple_frames(frame_server, imd_server, simple_frame_data, client):
    # wait for client to connect.
    time.sleep(0.1)
    frame_server.send_frame(0, simple_frame_data)
    frame_server.send_frame(1, simple_frame_data)
    time.sleep(0.5)
    assert len(client.frames) == 2


def test_reconnect_receive(frame_server, imd_server, client, simple_frame_data):
    frame_server.send_frame(0, simple_frame_data)
    client.close()
    assert client.latest_frame is None
    client.connect()
    time.sleep(0.2)
    assert client.latest_frame is not None


def test_close_interaction(frame_server, imd_server, interaction, client):
    id = client.start_interaction(interaction)
    time.sleep(0.2)
    assert len(imd_server.service.interactions) == 1
    client.close()
    time.sleep(0.2)
    assert len(imd_server.service.interactions) == 0


def test_running_imd(client):
    assert client.running_imd


def test_stop_interaction(frame_server, imd_server, client, interaction):
    id = client.start_interaction(interaction)
    time.sleep(0.2)
    assert len(imd_server.service.interactions) == 1
    client.stop_interaction(id)
    time.sleep(0.2)
    assert len(imd_server.service.interactions) == 0


def test_start_interaction(client, frame_server, imd_server, interaction):
    client.start_interaction(interaction)
    time.sleep(0.2)
    assert len(imd_server.service.interactions) == 1


def test_update_interaction(frame_server, imd_server, interaction, client):
    id = client.start_interaction(interaction)
    interaction.position = [2, 2, 2]
    client.update_interaction(id, interaction)
    time.sleep(0.5)
    assert len(imd_server.service.interactions) == 1
    assert np.allclose(list(imd_server.service.interactions.values())[0].position, (2, 2, 2))


def test_no_imd(frame_server, interaction):
    client = NarupaClient(run_imd=False)
    with pytest.raises(ValueError):
        client.start_interaction(interaction)
    client.close()


def test_no_imd_update(frame_server, interaction):
    client = NarupaClient(run_imd=False)
    with pytest.raises(ValueError):
        client.update_interaction(0, interaction)
    client.close()


def test_no_imd_stop(frame_server, interaction):
    client = NarupaClient(run_imd=False)
    with pytest.raises(ValueError):
        client.stop_interaction(0)
    client.close()
