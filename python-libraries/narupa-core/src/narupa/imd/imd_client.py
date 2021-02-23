# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.
import logging
from uuid import uuid4
from typing import Dict, Set

import grpc
from narupa.core import NarupaClient
from narupa.imd.particle_interaction import ParticleInteraction
from narupa.imd.imd_state import (
    INTERACTION_PREFIX,
    interaction_to_dict,
    dict_to_interaction,
)
from narupa.utilities.change_buffers import DictionaryChange


class ImdClient(NarupaClient):
    """
    A simple IMD client, primarily for testing the IMD server.

    """
    _local_interaction_ids: Set[str]
    _logger: logging.Logger

    def __init__(self, *,
                 channel: grpc.Channel,
                 make_channel_owner: bool = False):
        super().__init__(channel=channel, make_channel_owner=make_channel_owner)
        self._local_interaction_ids = set()
        self._logger = logging.getLogger(__name__)

    @property
    def interactions(self) -> Dict[str, ParticleInteraction]:
        with self.lock_state() as state:
            return {
                key: dict_to_interaction(value)
                for key, value in state.items()
                if key.startswith(INTERACTION_PREFIX)
            }

    def start_interaction(self) -> str:
        """
        Start an interaction

        :return: A unique identifier for the interaction.
        """
        interaction_id = INTERACTION_PREFIX + str(uuid4())
        self._local_interaction_ids.add(interaction_id)
        return interaction_id

    def update_interaction(
            self,
            interaction_id: str,
            interaction: ParticleInteraction,
    ):
        """
        Updates the interaction identified with the given interaction_id on the server with
        parameters from the given interaction.

        :param interaction_id: The unique interaction ID, created with
            :func:`~ImdClient.start_interaction`, that identifies the
            interaction to update.
        :param interaction: The :class: ParticleInteraction providing new
            parameters for the interaction.

        :raises: KeyError, if the given interaction ID does not exist.
        """
        if interaction_id not in self._local_interaction_ids:
            raise KeyError("Attempted to update an interaction with an "
                           "unknown interaction ID.")
        change = DictionaryChange(updates={
            interaction_id: interaction_to_dict(interaction),
        })
        return self.attempt_update_state(change)

    def stop_interaction(self, interaction_id: str) -> bool:
        """
        Stops the interaction identified with the given interaction_id on the server.

        :param interaction_id: The unique interaction ID, created with
            :func:`~ImdClient.start_interaction`, that identifies the
            interaction to stop.

        :raises: KeyError, if the given interaction ID does not exist.
        """
        if interaction_id not in self._local_interaction_ids:
            raise KeyError("Attempted to stop an interaction with an unknown "
                           "interaction ID.")
        self._local_interaction_ids.remove(interaction_id)
        change = DictionaryChange(
            removals=set([interaction_id])
        )
        return self.attempt_update_state(change)

    def stop_all_interactions(self):
        """
        Stops all active interactions governed by this client.
        """
        for interaction_id in list(self._local_interaction_ids):
            self.stop_interaction(interaction_id)

    def close(self):
        """
        Closes the IMD client.
        """
        try:
            self.stop_all_interactions()
        except grpc.RpcError as e:
            self._logger.exception(e)
            raise e
        finally:
            super().close()
