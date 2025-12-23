import numpy
import pytest
from nanover.app.multiuser import (
    RADIAL_ORIENT_COMMAND_KEY,
    MULTIUSER_ORIGIN_PREFIX,
)
from nanover.testing import assert_in_soon
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

    user_ids = ["test" + str(i) for i in range(avatar_count)]

    update = DictionaryChange({"avatar." + id: [] for id in user_ids})
    client.update_state(update)

    # wait for state update
    if user_ids:
        assert_in_soon(
            lambda: next(iter(update.updates)),
            lambda: server.state_dictionary.copy_content(),
        )

    client.run_command_blocking(RADIAL_ORIENT_COMMAND_KEY, radius=radius)

    # wait for state update
    if user_ids:
        assert_true_soon(
            lambda: [
                key
                for key in client.latest_multiplayer_values
                if key.startswith(MULTIUSER_ORIGIN_PREFIX)
            ],
        )

    origin_keys = {MULTIUSER_ORIGIN_PREFIX + id for id in user_ids}
    assert origin_keys <= client.latest_multiplayer_values.keys()

    positions = [
        value["position"]
        for key, value in client.latest_multiplayer_values.items()
        if key.startswith(MULTIUSER_ORIGIN_PREFIX)
    ]

    assert all(numpy.linalg.norm(position) == radius for position in positions)
