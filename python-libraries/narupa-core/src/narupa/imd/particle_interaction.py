# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.
"""
Module providing a wrapper class around the protobuf interaction message.
"""
from numbers import Number
from typing import Collection

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
    _interaction: imd_pb2.ParticleInteraction

    TYPE_KEY = "type"
    SCALE_KEY = "scale"
    MASS_WEIGHTED_KEY = "mass_weighted"
    RESET_VELOCITIES_KEY = "reset_velocities"
    MAX_FORCE_KEY = "max_force"

    def __init__(self,
                 player_id: str,
                 interaction_id: str,
                 position=(0, 0, 0),
                 particles=(),
                 interaction_type='gaussian',
                 scale=1,
                 mass_weighted=True,
                 reset_velocities=False,
                 max_force=20000):
        self._interaction = imd_pb2.ParticleInteraction(player_id=player_id, interaction_id=interaction_id)
        self.position = position
        self._properties = self._interaction.properties
        self.scale = scale
        self.type = interaction_type
        self.mass_weighted = mass_weighted
        self.reset_velocities = reset_velocities
        self.particles = particles
        self.max_force = max_force

    @classmethod
    def from_proto(cls,
                   interaction_proto: imd_pb2.ParticleInteraction,
                   default_interaction_type='gaussian',
                   default_scale=1,
                   default_mass_weighted=True,
                   default_reset_velocities=False,
                   default_max_force=20000):
        """
        Initialises an interaction from the protobuf representation.

        :param interaction_proto: The protobuf representation of the interaction.
        """
        interaction = cls(
            player_id=interaction_proto.player_id,
            interaction_id=interaction_proto.interaction_id,
        )
        interaction._interaction = interaction_proto
        interaction._properties = interaction_proto.properties
        set_default_property(interaction.properties, cls.TYPE_KEY, default_interaction_type)
        set_default_property(interaction.properties, cls.MASS_WEIGHTED_KEY, default_mass_weighted)
        set_default_property(interaction.properties, cls.SCALE_KEY, default_scale)
        set_default_property(interaction.properties, cls.RESET_VELOCITIES_KEY, default_reset_velocities)
        set_default_property(interaction.properties, cls.MAX_FORCE_KEY, default_max_force)

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

        :return: The player ID associated with this interaction.
        """
        return self._interaction.player_id

    @property
    def interaction_id(self) -> str:
        """
        Gets the interaction ID associated with this interaction.

        :return: The interaction ID associated with this interaction.
        """
        return self._interaction.interaction_id

    @property
    def type(self) -> str:
        """
        Gets the type of interaction being applied, default 'gaussian'.

        :return: The type of interaction being applied.
        """
        return self._get_property('type')

    @type.setter
    def type(self, value: str):
        """
        Sets the interaction type.

        :param value: Interaction type to apply. Typically 'gaussian' or 'spring'.
        """
        self._set_property('type', value)

    @property
    def scale(self) -> float:
        """
        Gets the scale factor of the interaction, which defaults to 1.

        Adjusting this changes the strength of the interactive force applied.

        :return: The scale factor of the interaction.
        """
        return self._get_property('scale')

    @scale.setter
    def scale(self, value: Number):
        """
        Sets the scale factor of the interaction.

        :param value: The new scale factor to set.
        """
        self._set_property('scale', value)

    @property
    def position(self) -> np.array:
        """
        Gets the position of the interaction, which defaults to ``[0,0,0]``

        :return: The position of the interaction, in nanometers.
        """
        return np.array(self._interaction.position)

    @position.setter
    def position(self, position: Collection[float]):
        """
        Set the position of the interaction

        :param position: 3 dimensional vector position of interaction, in nanometers.
        """
        if len(position) != 3:
            raise ValueError(f"Tried to set an interaction position that did not have 3 points. "
                             f"It had {len(position)} points.")
        self._interaction.position[:] = position

    @property
    def particles(self) -> np.ndarray:
        """
        Gets the list of particles this interaction applies to.

        :return: The list of the indices of the particles this interaction applies to.
        """
        return np.array(self._interaction.particles)

    @particles.setter
    def particles(self, particles: Collection[int]):
        """
        Set the particles of the interaction.

        :param particles: A collection of particles. If it contains duplicates, these will be removed.
        """
        self._interaction.particles[:] = np.unique(particles)

    @property
    def max_force(self) -> float:
        """
        Gets the maximum force this interaction will be allowed to apply to the system.

        :return: The maximum energy, in kJ/(mol*nm), the interaction will be allowed to apply to the system.
        """
        return self._get_property(self.MAX_FORCE_KEY)

    @max_force.setter
    def max_force(self, value: float):
        """
        Sets the maximum force this interaction will be allowed to apply to the system.

        :param value: New maximum force, in kJ/(mol*nm).
        """
        self._set_property(self.MAX_FORCE_KEY, value)

    @property
    def mass_weighted(self) -> bool:
        """
        Indicates whether this interaction should be mass weighted, default `True`.

        :return: Whether to mass weight this interaction.
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

        :param value: Boolean value to set.
        """
        self._properties['mass_weighted'] = value

    @property
    def reset_velocities(self) -> bool:
        """
        Indicates whether this interaction should be reset the velocities of
        the atoms it interacts with after interaction, defaulting to False.

        :return: Whether to reset velocities after this interaction.
        """
        # TODO should we update these to set the property if it does not exist, for serialisation?
        try:
            result = self._properties[self.RESET_VELOCITIES_KEY]
        except ValueError:
            return False
        else:
            return result

    @reset_velocities.setter
    def reset_velocities(self, value):
        """
        Sets this interaction to reset the velocities of the selected atoms or
        not after the interaction is complete.

        :param value: Whether to reset velocities after this interaction.
        """
        self._properties[self.RESET_VELOCITIES_KEY] = value

    @property
    def properties(self) -> Struct:
        """
        Gets the properties Struct field of the interaction structure.
        """
        return self._properties

    def _get_property(self, key):
        return self.properties[key]

    def _set_property(self, key, value):
        self._properties[key] = value
