from typing import List, Collection

import numpy as np

from narupa.protocol.imd.imd_pb2 import Interaction as InteractionGrpc


class Interaction:
    """
    A wrapper around the protobuf representation of an interaction.
    Provides easy to use getters and setters.

    For convenience, the getters all copy the underlying data into numpy arrays,
    rather than the low level containers used by protobuf.

    :param player_id: The player ID to be associated with the interaction.

    """
    _interaction: InteractionGrpc

    def __init__(self, player_id: str = "1"):
        """

        """
        self._interaction = InteractionGrpc(player_id=player_id)
        self._interaction.position[:] = [0, 0, 0]

    @property
    def player_id(self) -> str:
        """
        Gets the player ID associated with this interaction.
        :return:
        """
        return self._interaction.player_id

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
    def atoms(self) -> Collection:
        """
        Gets the list of atoms this interaction applies to.
        :return:
        """
        return np.array(self._interaction.atoms)

    @atoms.setter
    def atoms(self, atoms: Collection):
        """
        Set the atoms of the interaction.
        :param atoms: A set of atoms. If it contains duplicates, these will be removed.
        :return:
        """
        self._interaction.atoms[:] = np.unique(atoms)
