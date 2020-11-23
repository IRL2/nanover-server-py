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
        updates = {
            MULTIUSER_ORIGIN_PREFIX + avatar_id: {
                "position": [0, 0, 0],
                "rotation": [0, 0, 0, 1],
            }
            for avatar_id in avatar_ids
        }
        server.update_state(None, DictionaryChange(updates))

    server.register_command(RADIAL_ORIENT_COMMAND_KEY, radially_orient)
