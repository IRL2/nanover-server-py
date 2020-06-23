import time

import pytest
from narupa.app import NarupaImdClient
from narupa.core import NarupaServer

UPDATE_TIME = 0.01


def assert_selections_base_equal(selection, name, particle_ids):
    """
    Assert that a selection as the expected name and selected particles.

    The method test both the name and the selected particles in one go.
    """
    selection_tuple = (selection.selection_name, selection.selected_particle_ids)
    expectation_tuple = (name, particle_ids)
    assert selection_tuple == expectation_tuple


def assert_number_and_get_first_selection(client, expected_number_of_selections):
    """
    Assert that a client returns the expected number of selections and return
    the first one.

    This function can only be used if the expected number of selections is a
    strictly positive.
    """
    if expected_number_of_selections <= 0:
        raise ValueError('The expected number of selections must be strictly positive.')
    selections = list(client.selections)
    assert len(selections) == expected_number_of_selections
    return selections[0]


@pytest.fixture
def server_clients():
    """
    Provides a multiplayer server hosting on an available port on localhost,
    and two Narupa clients connected to it.
    """
    server = NarupaServer(address='localhost', port=0)
    client1 = NarupaImdClient(multiplayer_address=(server.address, server.port))
    client2 = NarupaImdClient(multiplayer_address=(server.address, server.port))

    client1.subscribe_multiplayer(interval=0)
    client2.subscribe_multiplayer(interval=0)

    with server, client1, client2:
        yield server, client1, client2


@pytest.fixture
def server_clients_with_selection(server_clients):
    server, client1, client2 = server_clients
    client1.create_selection('Selection 1', [0, 1, 2])
    time.sleep(UPDATE_TIME)
    yield server, client1, client2


def test_client_key_sharing(server_clients):
    server, client1, client2 = server_clients

    value = 'value'
    client1.set_shared_value('key', value)
    time.sleep(UPDATE_TIME)
    assert client1.get_shared_value('key') == value


def test_create_selection(server_clients):
    server, client1, client2 = server_clients

    client1.create_selection('Selection 1', [0, 1, 2])

    time.sleep(UPDATE_TIME)

    selection = assert_number_and_get_first_selection(client2, 1)
    assert_selections_base_equal(selection, 'Selection 1', {0, 1, 2})


def test_create_empty_selection(server_clients):
    server, client1, client2 = server_clients

    client1.create_selection('Empty Selection')

    time.sleep(UPDATE_TIME)
    selection = assert_number_and_get_first_selection(client2, 1)
    assert_selections_base_equal(selection, 'Empty Selection', set())


def test_update_selection(server_clients_with_selection):
    server, client1, client2 = server_clients_with_selection

    selection = assert_number_and_get_first_selection(client2, 1)
    selection.selection_name = "Selection 2"
    selection.selected_particle_ids = {3, 4, 5}
    client2.update_selection(selection)

    time.sleep(UPDATE_TIME)
    selection = assert_number_and_get_first_selection(client1, 1)
    assert_selections_base_equal(selection, 'Selection 2', {3, 4, 5})


def test_update_selection_with(server_clients_with_selection):
    server, client1, client2 = server_clients_with_selection

    selection = assert_number_and_get_first_selection(client2, 1)
    with selection.modify():
        selection.selection_name = "Selection 2"
        selection.selected_particle_ids = {3, 4, 5}

    time.sleep(UPDATE_TIME)
    assert_number_and_get_first_selection(client1, 1)
    assert_selections_base_equal(selection, 'Selection 2', {3, 4, 5})


def test_remove_selection(server_clients_with_selection):
    server, client1, client2 = server_clients_with_selection

    selection = assert_number_and_get_first_selection(client2, 1)
    client2.remove_selection(selection)

    time.sleep(UPDATE_TIME)
    assert len(list(client1.selections)) == 0


def test_remove_selection_remove_method(server_clients_with_selection):
    server, client1, client2 = server_clients_with_selection

    selection = assert_number_and_get_first_selection(client2, 1)
    selection.remove()

    time.sleep(UPDATE_TIME)
    assert len(list(client1.selections)) == 0


def test_clear_selections(server_clients_with_selection):
    server, client1, client2 = server_clients_with_selection

    client1.create_selection('Selection 2', [3, 4])

    time.sleep(UPDATE_TIME)

    assert_number_and_get_first_selection(client1, 2)
    assert_number_and_get_first_selection(client2, 2)

    client1.clear_selections()
    time.sleep(UPDATE_TIME)

    assert len(list(client1.selections)) == 0
    assert len(list(client2.selections)) == 0


def test_get_selection(server_clients_with_selection):
    server, client1, client2 = server_clients_with_selection

    selection = assert_number_and_get_first_selection(client1, 1)
    id_ = selection.selection_id

    time.sleep(UPDATE_TIME)

    selection = client2.get_selection(id_)
    assert_selections_base_equal(selection, 'Selection 1', {0, 1, 2})


def test_get_selection_missing(server_clients_with_selection):
    server, client1, client2 = server_clients_with_selection
    with pytest.raises(KeyError):
        _ = client2.get_selection("selection.invalid_id")


def test_root_selection_fields(server_clients):
    server, client1, client2 = server_clients

    selection = client1.root_selection
    assert not selection.selected_particle_ids


def test_root_selection_set_field(server_clients):
    server, client1, client2 = server_clients

    selection = client1.root_selection
    with selection.modify():
        selection.hide = True

    time.sleep(UPDATE_TIME)
    selection = client2.root_selection
    assert selection.hide is True


def test_remove_selection_while_in_use(server_clients_with_selection):
    server, client1, client2 = server_clients_with_selection

    selection = assert_number_and_get_first_selection(client1, 1)
    with selection.modify():
        selection.selection_name = "Selection 2"
        selection_from_2 = assert_number_and_get_first_selection(client2, 1)
        assert_selections_base_equal(selection_from_2, 'Selection 1', {0, 1, 2})

        client2.clear_selections()
        time.sleep(UPDATE_TIME)

        assert len(list(client1.selections)) == 0

    time.sleep(UPDATE_TIME)

    assert len(list(client1.selections)) == 1
    assert len(list(client2.selections)) == 1
