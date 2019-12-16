import time

import pytest
from google.protobuf.struct_pb2 import Value
from narupa.app import NarupaImdClient
from narupa.multiplayer import MultiplayerServer

UPDATE_TIME = 0.05


@pytest.fixture
def server_clients():
    """
    Provides a multiplayer server hosting on an available port on localhost,
    and two Narupa clients connected to it.
    """
    server = MultiplayerServer(address='localhost', port=0)
    client1 = NarupaImdClient(
        multiplayer_port=server.port,
        address=server.address
    )
    client2 = NarupaImdClient(
        multiplayer_port=server.port,
        address=server.address
    )

    client1.join_multiplayer("player 1")
    client2.join_multiplayer("player 2")

    time.sleep(0.02)

    with server, client1, client2:
        yield server, client1, client2


def test_client_key_sharing(server_clients):
    server, client1, client2 = server_clients

    value = Value(string_value='value')
    client1.set_shared_value('key', value)
    time.sleep(UPDATE_TIME)
    assert client1.get_shared_value('key') == value


def test_create_selection(server_clients):
    server, client1, client2 = server_clients

    client1.create_selection('Selection 1', [0, 1, 2])

    time.sleep(UPDATE_TIME)

    assert len(list(client2.selections)) == 1
    selection = list(client2.selections)[0]
    assert selection.selection_name == 'Selection 1'
    assert selection.selected_particle_ids == {0, 1, 2}


def test_create_empty_selection(server_clients):
    server, client1, client2 = server_clients

    client1.create_selection('Empty Selection')

    time.sleep(UPDATE_TIME)

    assert len(list(client2.selections)) == 1
    selection = list(client2.selections)[0]
    assert selection.selected_particle_ids == set()


def test_update_selection(server_clients):
    server, client1, client2 = server_clients

    client1.create_selection('Selection 1', [0, 1, 2])

    time.sleep(UPDATE_TIME)

    assert len(list(client2.selections)) == 1
    selection = list(client2.selections)[0]
    assert selection.selection_name == 'Selection 1'
    assert selection.selected_particle_ids == {0, 1, 2}

    selection.selection_name = "Selection 2"
    selection.selected_particle_ids = {3, 4, 5}
    client2.update_selection(selection)

    time.sleep(UPDATE_TIME)

    assert len(list(client1.selections)) == 1
    selection = list(client1.selections)[0]
    assert selection.selection_name == 'Selection 2'
    assert selection.selected_particle_ids == {3, 4, 5}


def test_update_selection_with(server_clients):
    server, client1, client2 = server_clients

    client1.create_selection('Selection 1', [0, 1, 2])

    time.sleep(UPDATE_TIME)

    assert len(list(client2.selections)) == 1
    selection = list(client2.selections)[0]
    assert selection.selection_name == 'Selection 1'
    assert selection.selected_particle_ids == {0, 1, 2}

    with selection.modify():
        selection.selection_name = "Selection 2"
        selection.selected_particle_ids = {3, 4, 5}

    time.sleep(UPDATE_TIME)

    assert len(list(client1.selections)) == 1
    selection = list(client1.selections)[0]
    assert selection.selection_name == 'Selection 2'
    assert selection.selected_particle_ids == {3, 4, 5}


def test_remove_selection(server_clients):
    server, client1, client2 = server_clients

    client1.create_selection('Selection 1', [0, 1, 2])

    time.sleep(UPDATE_TIME)

    assert len(list(client2.selections)) == 1

    selection = list(client2.selections)[0]
    client2.remove_selection(selection)

    time.sleep(UPDATE_TIME)

    assert len(list(client1.selections)) == 0


def test_remove_selection_remove_method(server_clients):
    server, client1, client2 = server_clients

    client1.create_selection('Selection 1', [0, 1, 2])

    time.sleep(UPDATE_TIME)

    assert len(list(client2.selections)) == 1

    selection = list(client2.selections)[0]
    selection.remove()

    time.sleep(UPDATE_TIME)

    assert len(list(client1.selections)) == 0


def test_clear_selections(server_clients):
    server, client1, client2 = server_clients

    client1.create_selection('Selection 1', [0, 1, 2])

    time.sleep(UPDATE_TIME)

    client1.create_selection('Selection 2', [3, 4])

    time.sleep(UPDATE_TIME)

    assert len(list(client1.selections)) == 2
    assert len(list(client2.selections)) == 2

    client1.clear_selections()

    time.sleep(UPDATE_TIME)

    assert len(list(client1.selections)) == 0
    assert len(list(client2.selections)) == 0


def test_get_selection(server_clients):
    server, client1, client2 = server_clients

    selection = client1.create_selection('Selection 1', [0, 1, 2])
    id = selection.selection_id

    time.sleep(UPDATE_TIME)

    assert len(list(client2.selections)) == 1

    assert client2.get_selection(id) is not None

    with pytest.raises(KeyError):
        _ = client2.get_selection("selection.invalid_id")


def test_root_selection_fields(server_clients):
    server, client1, client2 = server_clients

    selection = client1.root_selection

    assert selection is not None

    assert selection.selected_particle_ids is None


def test_root_selection_set_field(server_clients):
    server, client1, client2 = server_clients

    selection = client1.root_selection

    assert selection is not None

    with selection.modify():
        selection.hide = True

    time.sleep(UPDATE_TIME)

    selection = client2.root_selection

    assert selection.hide is True


def test_remove_selection_while_in_use(server_clients):
    server, client1, client2 = server_clients

    selection = client1.create_selection("Selection 1")

    assert selection is not None

    time.sleep(UPDATE_TIME)

    with selection.modify():
        selection.selection_name = "Selection 2"

        assert len(list(client2.selections)) == 1

        client2.clear_selections()

        time.sleep(UPDATE_TIME)

        assert len(list(client1.selections)) == 0

    time.sleep(UPDATE_TIME)

    assert len(list(client1.selections)) == 1
    assert len(list(client2.selections)) == 1
