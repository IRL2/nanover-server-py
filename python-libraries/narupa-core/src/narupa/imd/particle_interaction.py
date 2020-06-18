# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.
"""
Module providing a wrapper class around the protobuf interaction message.
"""
import math
from typing import Collection, Dict, Any, Iterable

import numpy as np
from google.protobuf.struct_pb2 import Struct

import narupa.protocol.imd.imd_pb2 as imd_pb2
from narupa.utilities.protobuf_utilities import dict_to_struct, struct_to_dict

DEFAULT_MAX_FORCE = 20000


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
    :param interaction_id: The interaction ID to be associated with the
        interaction. Typically, this identifies the VR controller, or other
        input device.
    :param interaction_type: The type of interaction being used, default is
        'gaussian' for a Gaussian force.
    :param scale: The scale factor applied to the interaction, default is 1.
    :param mass_weighted: Whether the interaction will be mass weighted or not.
    :param reset_velocities: Whether to reset velocities after interacting.
    :param max_force: The maximum force that will be allowed to be applied to a given atom in a given cartesian
        direction. Helps maintain stability for unbounded potentials.

    """

    TYPE_KEY = "type"
    SCALE_KEY = "scale"
    MASS_WEIGHTED_KEY = "mass_weighted"
    RESET_VELOCITIES_KEY = "reset_velocities"
    MAX_FORCE_KEY = "max_force"

    def __init__(self,
                 player_id: str,
                 interaction_id: str,
                 position=(0., 0., 0.),
                 particles=(),
                 interaction_type='gaussian',
                 scale=1,
                 mass_weighted=True,
                 reset_velocities=False,
                 max_force=DEFAULT_MAX_FORCE,
                 **kwargs):
        self.player_id = player_id
        self.interaction_id = interaction_id
        self.position = position
        self.particles = particles
        self.scale = scale
        self.type = interaction_type
        self.mass_weighted = mass_weighted
        self.reset_velocities = reset_velocities
        self.max_force = max_force
        self.properties = dict(kwargs)

    @classmethod
    def from_proto(cls,
                   interaction_proto: imd_pb2.ParticleInteraction,
                   default_interaction_type='gaussian',
                   default_scale=1,
                   default_mass_weighted=True,
                   default_reset_velocities=False,
                   default_max_force=DEFAULT_MAX_FORCE):
        """
        Initialises an interaction from the protobuf representation.

        :param interaction_proto: The protobuf representation of the interaction.
        """

        properties_dict = struct_to_dict(interaction_proto.properties)

        type = properties_dict.pop(cls.TYPE_KEY, default_interaction_type)

        mass_weighted = default_mass_weighted
        if cls.MASS_WEIGHTED_KEY in interaction_proto.properties:
            mass_weighted = properties_dict[cls.MASS_WEIGHTED_KEY]
            del properties_dict[cls.MASS_WEIGHTED_KEY]

        scale = default_scale
        if cls.SCALE_KEY in interaction_proto.properties:
            scale = properties_dict[cls.SCALE_KEY]
            del properties_dict[cls.SCALE_KEY]

        reset_velocities = default_reset_velocities
        if cls.RESET_VELOCITIES_KEY in interaction_proto.properties:
            reset_velocities = properties_dict[cls.RESET_VELOCITIES_KEY]
            del properties_dict[cls.RESET_VELOCITIES_KEY]

        max_force = default_max_force
        if cls.MAX_FORCE_KEY in interaction_proto.properties:
            max_force = properties_dict[cls.MAX_FORCE_KEY]
            del properties_dict[cls.MAX_FORCE_KEY]

        interaction = cls(
            player_id=interaction_proto.player_id,
            interaction_id=interaction_proto.interaction_id,
            position=interaction_proto.position,
            particles=interaction_proto.particles,
            interaction_type=type,
            mass_weighted=mass_weighted,
            scale=scale,
            reset_velocities=reset_velocities,
            max_force=max_force,
            **properties_dict
        )

        return interaction

    @property
    def proto(self) -> imd_pb2.ParticleInteraction:
        """
        Gets the underlying protobuf representation.

        :return: The underlying protobuf Interaction representation.
        """
        interaction = imd_pb2.ParticleInteraction()
        interaction.position.extend(self.position)
        interaction.particles.extend(self.particles)
        interaction.player_id = self.player_id
        interaction.interaction_id = self.interaction_id
        properties = {}
        for key, value in self._properties.items():
            properties[key] = value
        properties[self.TYPE_KEY] = self.type
        properties[self.MASS_WEIGHTED_KEY] = self.mass_weighted
        properties[self.SCALE_KEY] = self.scale
        properties[self.RESET_VELOCITIES_KEY] = self.reset_velocities
        properties[self.MAX_FORCE_KEY] = self.max_force
        for key, value in dict_to_struct(properties).items():
            interaction.properties[key] = value
        return interaction

    @property
    def player_id(self) -> str:
        """
        The player ID associated with this interaction.
        """
        return self._player_id

    @player_id.setter
    def player_id(self, value: str):
        self._player_id = value

    @property
    def interaction_id(self) -> str:
        """
        The interaction ID associated with this interaction.
        """
        return self._interaction_id

    @interaction_id.setter
    def interaction_id(self, value: str):
        self._interaction_id = value

    @property
    def type(self) -> str:
        """
        The type of interaction being applied, default 'gaussian'.
        """
        return self._type

    @type.setter
    def type(self, value: str):
        self._type = value

    @property
    def scale(self) -> float:
        """
        The scale factor of the interaction, which defaults to 1.

        Adjusting this changes the strength of the interactive force applied.
        """
        return self._scale

    @scale.setter
    def scale(self, value: float):
        self._scale = value

    @property
    def position(self) -> np.array:
        """
        The position of the interaction in nanometers, which defaults to ``[0,0,0]``
        """
        return self._position

    @position.setter
    def position(self, position: Iterable[float]):
        converted = np.array(position)
        if len(converted) != 3:
            raise ValueError(f"Position expected 3d vector, instead received: {position}")
        self._position = converted

    @property
    def particles(self) -> np.ndarray:
        """
        The list of particles this interaction applies to.
        """
        return self._particles

    @particles.setter
    def particles(self, particles: Collection[int]):
        self._particles = np.unique(particles)

    @property
    def max_force(self) -> float:
        """
        The maximum force, in kJ/(mol*nm), this interaction will be allowed to apply to the system.
        """
        return self._max_force

    @max_force.setter
    def max_force(self, value: float):
        if math.isnan(value):
            value = math.inf
        self._max_force = value

    @property
    def mass_weighted(self) -> bool:
        """
        Indicates whether this interaction should be mass weighted, default `True`.
        """
        return self._mass_weighted

    @mass_weighted.setter
    def mass_weighted(self, value: bool):
        self._mass_weighted = value

    @property
    def reset_velocities(self) -> bool:
        """
        Indicates whether this interaction should be reset the velocities of
        the atoms it interacts with after interaction, defaulting to False.
        """
        return self._reset_velocities

    @reset_velocities.setter
    def reset_velocities(self, value: bool):
        self._reset_velocities = value

    @property
    def properties(self) -> Dict[str, Any]:
        """
        Gets the other properties for this interaction
        """
        return self._properties

    @properties.setter
    def properties(self, value: Dict[str, Any]):
        self._properties = value
