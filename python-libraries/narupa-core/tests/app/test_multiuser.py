import time
import numpy
import pytest
from narupa.app import NarupaApplicationServer, NarupaImdClient
from narupa.app.multiuser import (
    RADIAL_ORIENT_COMMAND_KEY, MULTIUSER_ORIGIN_PREFIX,
)
from narupa.utilities.change_buffers import DictionaryChange


@pytest.mark.parametrize('avatars', (0, 1, 4))
@pytest.mark.parametrize('radius', (0, 1, 5))
def test_radial_orient(avatars, radius):
    """
    Test that the radial orientation command creates the correct number of
    origin entries at the correct radius and spacing.
    """
    with NarupaApplicationServer.basic_server(port=0) as server:
        with NarupaImdClient.connect_to_single_server("localhost", server.port) as client:
            client.subscribe_multiplayer(interval=0)

            update = DictionaryChange({
                "avatar.test" + str(i): [] for i in range(avatars)
            })
            client.attempt_update_multiplayer_state(update)
            client.run_command(RADIAL_ORIENT_COMMAND_KEY, radius=radius)
            time.sleep(.001)

            origins = [
                value
                for key, value in client.latest_multiplayer_values.items()
                if key.startswith(MULTIUSER_ORIGIN_PREFIX)
            ]
            positions = [origin["position"] for origin in origins]

            assert len(origins) == avatars
            assert all(numpy.linalg.norm(position) == radius for position in positions)
