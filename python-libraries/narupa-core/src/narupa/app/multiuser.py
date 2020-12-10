# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.
"""
Module for addon server behaviour specific to the NarupaIMD multiuser
experience.
"""
import math
from narupa.core import NarupaServer
from narupa.utilities.change_buffers import DictionaryChange

RADIAL_ORIENT_COMMAND_KEY = 'multiuser/radially-orient-origins'
MULTIUSER_AVATAR_PREFIX = 'avatar.'
MULTIUSER_ORIGIN_PREFIX = 'user-origin.'


def add_multiuser_commands(server: NarupaServer):
    """
    Add server commands specific to the NarupaIMD multiuser experience.
    """
    def radially_orient(radius=1):
        """
        For each avatar present, add a suggested origin to the multiuser state
        for the corresponding user id. Distributes the origins in a circle
        around the true origin, each facing inwards, and each displaced radially
        by the given radius.
        """
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
        updates = {
            MULTIUSER_ORIGIN_PREFIX + avatar_id: {
                "position": [radius * math.cos(angle), 0, radius * math.sin(angle)],
                "rotation": _angle_axis_quaternion_y(-angle),
            }
            for avatar_id, angle in zip(avatar_ids, angles)
        }
        server.update_state(None, DictionaryChange(updates))

    server.register_command(RADIAL_ORIENT_COMMAND_KEY, radially_orient)


def _angle_axis_quaternion_y(angle):
    return [
        0,
        math.sin(angle * .5),
        0,
        math.cos(angle * .5),
    ]
