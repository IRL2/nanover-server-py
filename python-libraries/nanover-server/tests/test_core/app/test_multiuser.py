import numpy
import pytest
from nanover.app.multiuser import (
    RADIAL_ORIENT_COMMAND_KEY,
    MULTIUSER_ORIGIN_PREFIX,
)
from nanover.testing.asserts import assert_true_soon
from nanover.utilities.change_buffers import DictionaryChange

from ..app.test_client_selections import server_clients


@pytest.mark.parametrize("avatar_count", (0, 1, 4))
@pytest.mark.parametrize("radius", (0, 1, 5))
def test_radial_orient(server_clients, avatar_count, radius):
    """
    Test that the radial orientation command creates the correct number of
    origin entries at the correct radius and spacing.
    """
    server, client, _ = server_clients
    client.update_available_commands()

    user_ids = ["test" + str(i) for i in range(avatar_count)]

    update = DictionaryChange({"avatar." + id: [] for id in user_ids})
    client.update_state(update)
    client.run_command_blocking(RADIAL_ORIENT_COMMAND_KEY, radius=radius)

    origins = {
        key: value
        for key, value in client.latest_multiplayer_values.items()
        if key.startswith(MULTIUSER_ORIGIN_PREFIX)
    }
    positions = [origin["position"] for origin in origins.values()]

    assert_true_soon(
        lambda: all((MULTIUSER_ORIGIN_PREFIX + id) in origins for id in user_ids)
    )

    assert_true_soon(
        lambda: all(numpy.linalg.norm(position) == radius for position in positions)
    )
