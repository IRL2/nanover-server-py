from typing import Collection

import numpy as np

from narupa.core.utility.value import ValueMap
from narupa.protocol.imd.imd_pb2 import Interaction as InteractionGrpc


class Interaction:
    """
    A wrapper around the protobuf representation of an interaction.
    Provides easy to use getters and setters.

    For convenience, the getters all copy the underlying data into numpy arrays,
    rather than the low level containers used by protobuf.

    :param player_id: The player ID to be associated with the interaction.
    :param interaction_id: The interaction ID to be associated with the interaction. Typically, this identifies
    the VR controller, or other input device.

    """
    _interaction: InteractionGrpc

    def __init__(self, player_id: str = "1", interaction_id="2"):
        """

        """
        self._interaction = InteractionGrpc(player_id=player_id, interaction_id=interaction_id)
        self._interaction.position[:] = [0, 0, 0]
        self._properties = self._interaction.properties

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

    @property
    def properties(self):
        """
        Gets the properties field of the interaction structure.
        :return:
        """
        return self._properties

    @particles.setter
    def particles(self, particles: Collection):
        """
        Set the particles of the interaction.
        :param particles: A set of particles. If it contains duplicates, these will be removed.
        :return:
        """
        self._interaction.particles[:] = np.unique(particles)
