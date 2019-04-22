from numbers import Number
from typing import Collection

import numpy as np
from google.protobuf.struct_pb2 import Struct

import narupa.protocol.imd.imd_pb2 as imd_pb2


class Interaction:
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
    _interaction: imd_pb2.Interaction

    def __init__(self, player_id: str = "1", interaction_id="0", position=(0,0,0), interaction_type='gaussian', scale=1):
        """

        """
        self._interaction = imd_pb2.Interaction(player_id=player_id, interaction_id=interaction_id)
        self._interaction.position[:] = position
        self._properties = self._interaction.properties
        self.properties['scale'] = scale
        self.properties['type'] = interaction_type

    @classmethod
    def from_proto(cls, interaction_grpc):
        """
        Initialises an interaction from the protobuf representation.
        :param interaction_proto:
        :return:
        """
        interaction = cls()
        interaction._interaction = interaction_grpc
        interaction._properties = interaction_grpc.properties
        return interaction

    @property
    def proto(self) -> imd_pb2.Interaction:
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
    def type(self, value:str):
        """
        Sets the interaction type.
        :param value:
        :return:
        """
        self._set_property('type', value)

    @property
    def scale(self) -> Number:
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
    def particles(self) -> Collection:
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
