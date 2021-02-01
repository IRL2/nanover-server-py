import time
import numpy
import pytest
from narupa.app.multiuser import (
    RADIAL_ORIENT_COMMAND_KEY, MULTIUSER_ORIGIN_PREFIX, add_multiuser_commands,
)
from narupa.utilities.change_buffers import DictionaryChange

from ..app.test_client_selections import server_clients
from ..core.test_narupa_client_server_state import IMMEDIATE_REPLY_WAIT_TIME


@pytest.mark.parametrize('avatar_count', (0, 1, 4))
@pytest.mark.parametrize('radius', (0, 1, 5))
def test_radial_orient(server_clients, avatar_count, radius):
    """
    Test that the radial orientation command creates the correct number of
    origin entries at the correct radius and spacing.
    """
    server, client, _ = server_clients
    add_multiuser_commands(server)
    client.update_available_commands()

    user_ids = ["test" + str(i) for i in range(avatar_count)]

    update = DictionaryChange({
        "avatar." + id: [] for id in user_ids
    })
    client.attempt_update_multiplayer_state(update)
    client.run_command(RADIAL_ORIENT_COMMAND_KEY, radius=radius)
    time.sleep(IMMEDIATE_REPLY_WAIT_TIME)

    origins = {
        key: value
        for key, value in client.latest_multiplayer_values.items()
        if key.startswith(MULTIUSER_ORIGIN_PREFIX)
    }
    positions = [origin["position"] for origin in origins.values()]

    assert all((MULTIUSER_ORIGIN_PREFIX + id) in origins for id in user_ids)
    assert all(numpy.linalg.norm(position) == radius for position in positions)
