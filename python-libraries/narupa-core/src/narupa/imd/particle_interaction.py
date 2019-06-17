# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.
"""
Module providing a wrapper class around the protobuf interaction message.
"""
from numbers import Number
from typing import Collection, List

import numpy as np
from google.protobuf.struct_pb2 import Struct

import narupa.protocol.imd.imd_pb2 as imd_pb2


def set_default_property(properties: Struct, key, default):
    if key not in properties:
        properties[key] = default


class ParticleInteraction:
    """
    A wrapper around the protobuf representation of an interaction.
    Provides easy to use getters and setters.

    For convenience, the getters all copy the underlying data into numpy arrays,
    rather than the low level containers used by protobuf.

    :param player_id: The player ID to be associated with the interaction.
    :param interaction_id: The interaction ID to be associated with the interaction. Typically, this identifies
        the VR controller, or other input device.
    :param interaction_type: The type of interaction being used, default is Gaussian force.
    :param scale: The scale factor applied to the interaction, default is 1.

    """
    _interaction: imd_pb2.ParticleInteraction

    TYPE_KEY = "type"
    SCALE_KEY = "scale"
    MASS_WEIGHTED_KEY = "mass_weighted"


    def __init__(self, player_id: str = "1",
                 interaction_id="0",
                 position=(0, 0, 0),
                 particles=(),
                 interaction_type='gaussian',
                 scale=1,
                 mass_weighted=True):
        self._interaction = imd_pb2.ParticleInteraction(player_id=player_id, interaction_id=interaction_id)
        self.position = position
        self._properties = self._interaction.properties
        self.scale = scale
        self.type = interaction_type
        self.mass_weighted = mass_weighted
        self.particles = particles

    @classmethod
    def from_proto(cls, interaction_proto, default_interaction_type='gaussian', default_scale=1, default_mass_weighted=True):
        """
        Initialises an interaction from the protobuf representation.
        :param interaction_proto: The protobuf representation of the interaction.
        :return:
        """
        interaction = cls()
        interaction._interaction = interaction_proto
        interaction._properties = interaction_proto.properties
        set_default_property(interaction.properties, cls.TYPE_KEY, default_interaction_type)
        set_default_property(interaction.properties, cls.MASS_WEIGHTED_KEY, default_scale)
        set_default_property(interaction.properties, cls.SCALE_KEY, default_mass_weighted)

        return interaction

    @property
    def proto(self) -> imd_pb2.ParticleInteraction:
        """
        Gets the underlying protobuf representation.
        :return: The underlying protobuf Interaction representation.
        """
        return self._interaction

    @property
    def player_id(self) -> str:
        """
        Gets the player ID associated with this interaction.
        :return:
        """
        return self._interaction.player_id

    @property
    def interaction_id(self) -> str:
        """
        Gets the interaction ID associated with this interaction.
        :return:
        """
        return self._interaction.interaction_id

    @property
    def type(self) -> str:
        """
        Gets the type of interaction being applied, default 'gaussian'.
        :return:
        """
        return self._get_property('type')

    @type.setter
    def type(self, value: str):
        """
        Sets the interaction type.
        :param value:
        :return:
        """
        self._set_property('type', value)

    @property
    def scale(self) -> float:
        """
        Gets the scale factor of the interaction, which defaults to 1.
        :return:
        """
        return self._get_property('scale')

    @scale.setter
    def scale(self, value: Number):
        self._set_property('scale', value)

    @property
    def position(self) -> Collection:
        """
        Gets the position of the interaction, which defaults to [0,0,0]
        :return:
        """
        return np.array(self._interaction.position)

    @position.setter
    def position(self, position: Collection):
        """
        Set the position of the interaction
        :param position: Vector3 position of interaction.
        :return:
        """
        if len(position) != 3:
            raise ValueError(f"Tried to set an interaction position that did not have 3 points. "
                             f"It had {len(position)} points.")
        self._interaction.position[:] = position

    @property
    def particles(self) -> np.ndarray:
        """
        Gets the list of particles this interaction applies to.
        :return:
        """
        return np.array(self._interaction.particles)

    @particles.setter
    def particles(self, particles: Collection):
        """
        Set the particles of the interaction.
        :param particles: A set of particles. If it contains duplicates, these will be removed.
        :return:
        """
        self._interaction.particles[:] = np.unique(particles)

    @property
    def mass_weighted(self) -> bool:
        """
        Indicates whether this interaction should be mass weighted, default True.
        :return:
        """
        try:
            result = self._properties['mass_weighted']
        except ValueError:
            return True
        else:
            return result

    @mass_weighted.setter
    def mass_weighted(self, value: bool):
        """
        Sets this interaction to be mass weighted or not.
        :param value:
        :return:
        """
        self._properties['mass_weighted'] = value


    @property
    def properties(self) -> Struct:
        """
        Gets the properties Struct field of the interaction structure.
        :return:
        """
        return self._properties

    def _get_property(self, key):
        return self.properties[key]

    def _set_property(self, key, value):
        self._properties[key] = value
