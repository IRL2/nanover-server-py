# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.
"""
Module providing a wrapper class around the protobuf interaction message.
"""
import math
from typing import Dict, Any, Iterable, Collection
import numpy as np

DEFAULT_MAX_FORCE = 20000.0
DEFAULT_FORCE_TYPE = "gaussian"


class ParticleInteraction:
    """
    A wrapper around the protobuf representation of an interaction.
    Provides easy to use getters and setters.

    For convenience, the getters all copy the underlying data into numpy arrays,
    rather than the low level containers used by protobuf.

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
                 position=(0., 0., 0.),
                 particles=(),
                 interaction_type=DEFAULT_FORCE_TYPE,
                 scale=1.0,
                 mass_weighted=True,
                 reset_velocities=False,
                 max_force=DEFAULT_MAX_FORCE,
                 **kwargs):
        self.position = position
        self.particles = particles
        self.scale = scale
        self.interaction_type = interaction_type
        self.mass_weighted = mass_weighted
        self.reset_velocities = reset_velocities
        self.max_force = max_force
        self.properties = dict(kwargs)

    @property
    def interaction_type(self) -> str:
        """
        The type of interaction being applied, default 'gaussian'.
        """
        return self._type

    @interaction_type.setter
    def interaction_type(self, value: str):
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
    def particles(self, particles: Collection[int]):
        if len(particles) < 2:
            self._particles = np.array(particles)
        else:
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
                and self.reset_velocities == other.reset_velocities and self.interaction_type == other.interaction_type
                and self.properties == other.properties
        )

    def __repr__(self):
        return (
            f"<ParticleInteraction"
            f" position:{self.position}"
            f" particles:{self.particles}"
            f" reset_velocities:{self.reset_velocities}"
            f" scale:{self.scale}"
            f" mass_weighted:{self.mass_weighted}"
            f" max_force:{self.max_force}"
            f" type:{self.interaction_type}"
            f" other:{self.properties}"
            ">"
        )
