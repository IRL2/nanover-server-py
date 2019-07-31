# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.

"""
Module with helper functions for editing the multiplayer scene.
"""

from google.protobuf.struct_pb2 import Struct

#TODO replace this with a real wrapper class around the scene properties.
def get_default_scene_properties():
    """
    Initialises the default scene properties of the scene in the centre of the space, with no rotation and scale of 1.
    :return: Scene properties.
    """
    position = [0, 0, 0]
    rotation = [0, 0, 0, 1]
    scale = 1

    properties = Struct()
    properties["position"] = position
    properties["rotation"] = rotation
    properties["scale"] = scale

    return properties


