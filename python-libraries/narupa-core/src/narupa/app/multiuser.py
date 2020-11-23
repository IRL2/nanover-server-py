import math
from narupa.core import NarupaServer
from narupa.utilities.change_buffers import DictionaryChange

RADIAL_ORIENT_COMMAND_KEY = 'multiuser/radially-orient-origins'
MULTIUSER_AVATAR_PREFIX = 'avatar.'
MULTIUSER_ORIGIN_PREFIX = 'user-origin.'


def add_multiuser_commands(server: NarupaServer):
    def radially_orient():
        # find relevant avatar ids
        state = server.copy_state()
        avatar_ids = [
            key[len(MULTIUSER_AVATAR_PREFIX):]
            for key in state
            if key.startswith(MULTIUSER_AVATAR_PREFIX)
        ]
        # generate an origin for each avatar
        count = len(avatar_ids)
        angles = [i * math.pi * 2 / count for i in range(count)]
        radius = 1
        updates = {
            MULTIUSER_ORIGIN_PREFIX + avatar_id: {
                "position": [radius * math.cos(angle), 0, radius * math.sin(angle)],
                "rotation": [0, math.sin(angle), 0, math.cos(angle)],
            }
            for avatar_id, angle in zip(avatar_ids, angles)
        }
        server.update_state(None, DictionaryChange(updates))

    server.register_command(RADIAL_ORIENT_COMMAND_KEY, radially_orient)
