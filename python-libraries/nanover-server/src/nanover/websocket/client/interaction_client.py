from typing import Mapping
from uuid import uuid4

from nanover.imd import ParticleInteraction
from nanover.imd.imd_state import (
    dict_to_interaction,
    INTERACTION_PREFIX,
    interaction_to_dict,
)
from nanover.utilities.change_buffers import DictionaryChange
from nanover.websocket.client.base_client import WebsocketClient


class InteractionClient(WebsocketClient):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._local_interaction_ids = set()

    @property
    def interactions(self) -> dict[str, ParticleInteraction]:
        with self._state_dictionary.lock_state() as state:
            return {
                key: dict_to_interaction(value)
                for key, value in state.items()
                if key.startswith(INTERACTION_PREFIX)
                # We can have a misformatted interactions in the shared state.
                # Here we ignore them silently.
                # TODO: log a warning when this happens.
                and isinstance(value, Mapping)
            }

    def start_interaction(self, interaction: ParticleInteraction | None = None) -> str:
        """
        Start an interaction with the IMD server.

        :param interaction: An optional :class: ParticleInteraction with which
            to begin.
        :return: The unique interaction ID of this interaction, which can be
            used to update the interaction with
            :func:`~NanoverClient.update_interaction`.

        :raises: ValueError, if the there is no IMD connection available.
        """
        interaction_id = INTERACTION_PREFIX + str(uuid4())  # type: ignore
        self._local_interaction_ids.add(interaction_id)
        if interaction is not None:
            self.update_interaction(interaction_id, interaction)
        return interaction_id

    def update_interaction(
        self,
        interaction_id: str,
        interaction: ParticleInteraction,
    ):
        """
        Updates the interaction identified with the given interaction_id on
        the server with parameters from the given interaction.

        :param interaction_id: The unique id of the interaction to be updated.
        :param interaction: The :class: ParticleInteraction providing new
            parameters for the interaction.

        :raises: ValueError, if the there is no IMD connection available, or
            if invalid parameters are passed to the server.
        :raises KeyError: if the given interaction ID does not exist.
        """
        if interaction_id not in self._local_interaction_ids:
            raise KeyError(
                "Attempted to update an interaction with an " "unknown interaction ID."
            )
        change = DictionaryChange(
            updates={
                interaction_id: interaction_to_dict(interaction),
            }
        )
        return self.update_state(change)

    def stop_interaction(self, interaction_id: str) -> bool:
        """
        Stops the interaction identified with the given interaction_id on the server.
        This method blocks until the server returns a reply indicating that the
        interaction has stopped.

        :param interaction_id: The unique interaction ID, created with
            :func:`~NanoverClient.start_interaction`, that identifies the
            interaction to stop.
        :return: An :class:`InteractionEndReply`, which is an empty message indicating
            successful termination of the interaction.

        :raises ValueError: if the there is no IMD connection available, or
            if invalid parameters are passed to the server.
        :raises KeyError: if the given interaction ID does not exist.
        """
        if interaction_id not in self._local_interaction_ids:
            raise KeyError(
                "Attempted to stop an interaction with an unknown " "interaction ID."
            )
        self._local_interaction_ids.remove(interaction_id)
        change = DictionaryChange(removals={interaction_id})
        return self.update_state(change)

    def stop_all_interactions(self):
        """
        Stops all active interactions governed by this client.
        """
        for interaction_id in list(self._local_interaction_ids):
            self.stop_interaction(interaction_id)

    def close(self):
        self.stop_all_interactions()
        super().close()
