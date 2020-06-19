# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.
"""
Module providing a wrapper class around the protobuf interaction message.
"""
import math
from typing import Dict, Any, Iterable

import numpy as np
import narupa.protocol.imd.imd_pb2 as imd_pb2
from narupa.utilities.protobuf_utilities import dict_to_struct, struct_to_dict

DEFAULT_MAX_FORCE = 20000.0
DEFAULT_FORCE_TYPE = "gaussian"


class ParticleInteraction:
    """
    A wrapper around the protobuf representation of an interaction.
    Provides easy to use getters and setters.

    For convenience, the getters all copy the underlying data into numpy arrays,
    rather than the low level containers used by protobuf.

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
                 interaction_id: str,
                 position=(0., 0., 0.),
                 particles=(),
                 interaction_type=DEFAULT_FORCE_TYPE,
                 scale=1.0,
                 mass_weighted=True,
                 reset_velocities=False,
                 max_force=DEFAULT_MAX_FORCE,
                 **kwargs):
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
    def from_proto(cls, interaction_proto: imd_pb2.ParticleInteraction):
        """
        Initialises an interaction from the protobuf representation.

        :param interaction_proto: The protobuf representation of the interaction.
        """
        proto_key_to_keyword = {
            cls.TYPE_KEY: 'interaction_type',
            cls.MASS_WEIGHTED_KEY: 'mass_weighted',
            cls.SCALE_KEY: 'scale',
            cls.RESET_VELOCITIES_KEY: 'reset_velocities',
            cls.MAX_FORCE_KEY: 'max_force',
        }

        properties = struct_to_dict(interaction_proto.properties)
        kwargs = {
            proto_key_to_keyword.get(proto_key, proto_key): value
            for proto_key, value in properties.items()
        }

        interaction = cls(
            interaction_id=interaction_proto.interaction_id,
            position=interaction_proto.position,
            particles=interaction_proto.particles,
            **kwargs,
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
        if not math.isfinite(value):
            raise ValueError("Scale must be finite")
        self._scale = float(value)

    @property
    def position(self) -> np.array:
        """
        The position of the interaction in nanometers, which defaults to ``[0 0 0]``
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
    def particles(self, particles: Iterable[int]):
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
            raise ValueError("Max force cannot be nan")
        self._max_force = float(value)

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

    def __eq__(self, other):
        return (
            isinstance(other, ParticleInteraction) and np.equal(self.particles, other.particles).all()
            and np.isclose(self.position, other.position).all() and math.isclose(self.max_force, other.max_force)
            and self.mass_weighted == other.mass_weighted and math.isclose(self.scale, other.scale)
            and self.reset_velocities == other.reset_velocities and self.type == other.type
            and self.properties == other.properties and self.interaction_id == other.interaction_id
        )

    def __repr__(self):
        return (
            f"<ParticleInteraction"
            f" interaction_id:{self.interaction_id}"
            f" position:{self.position}"
            f" particles:{self.particles}"
            f" reset_velocities:{self.reset_velocities}"
            f" scale:{self.scale}"
            f" mass_weighted:{self.mass_weighted}"
            f" max_force:{self.max_force}"
            f" type:{self.type}"
            f" other:{self.properties}"
            ">"
        )
