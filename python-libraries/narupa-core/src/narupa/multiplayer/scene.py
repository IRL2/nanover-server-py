# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.

"""
Module with helper functions for editing the multiplayer scene.
"""

import narupa.protocol.multiplayer.multiplayer_pb2 as multiplayer_proto
from narupa.multiplayer.multiplayer_lock import MultiplayerObjectLock

#TODO replace this with a real wrapper class around the scene properties.
def get_default_scene_properties():
    """
    Initialises the default scene properties of the scene in the centre of the space, with no rotation and scale of 1.
    :return: Scene properties.
    """
    position = [0, 0, 0]
    rotation = [0, 0, 0, 1]
    scale = 1
    properties = multiplayer_proto.SceneProperties(position=position, rotation=rotation, scale=scale)
    return properties


