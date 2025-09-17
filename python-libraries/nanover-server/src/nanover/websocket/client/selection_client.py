from typing import Iterable
from uuid import uuid4

from nanover.app import RenderingSelection
from nanover.app.client import SELECTION_ROOT_ID, SELECTION_ROOT_NAME
from nanover.utilities.change_buffers import DictionaryChange
from nanover.websocket.client.base_client import WebsocketClient


class SelectionClient(WebsocketClient):
    _next_selection_id = 0

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._player_id = str(uuid4())

    @property  # type: ignore
    def root_selection(self) -> RenderingSelection:
        """
        Get the root selection, creating it if it does not exist yet.

        :return: The selection representing the root selection of the system
        """
        try:
            root_selection = self.get_selection(SELECTION_ROOT_ID)
        except KeyError:
            root_selection = self._create_selection_from_id_and_name(
                SELECTION_ROOT_ID, SELECTION_ROOT_NAME
            )
        root_selection.selected_particle_ids = set()
        return root_selection

    def create_selection(
        self,
        name: str,
        particle_ids: Iterable[int] | None = None,
    ) -> RenderingSelection:
        """
        Create a particle selection with the given name.

        :param name: The user-friendly name of the selection.
        :param particle_ids: The indices of the particles to include in the selection.
        :return: The selection that was created.
        """
        if particle_ids is None:
            particle_ids = set()

        # Give the selection an ID based upon the multiplayer player ID and an incrementing counter
        selection_id = f"selection.{self._player_id}.{self._next_selection_id}"
        self._next_selection_id += 1

        # Create the selection and setup the particles that it contains
        selection = self._create_selection_from_id_and_name(selection_id, name)
        selection.set_particles(particle_ids)

        # Mark the selection as needing updating, which adds it to the shared value store.
        self.update_selection(selection)

        return selection

    def update_selection(self, selection: RenderingSelection):
        """
        Applies changes to the given selection to the shared key store.

        :param selection: The selection to update.
        """
        self.update_state(
            DictionaryChange(
                updates={selection.selection_id: selection.to_dictionary()}
            )
        )

    def remove_selection(self, selection: RenderingSelection):
        """
        Delete the given selection
        """
        self.update_state(DictionaryChange(removals={selection.selection_id}))

    def clear_selections(self):
        """
        Remove all selections in the system
        """
        selections = list(self.selections)
        for selection in selections:
            self.remove_selection(selection)

    @property  # type: ignore
    def selections(self) -> Iterable[RenderingSelection]:
        """
        Get all selections which are stored in the shared key store.

        :return: An iterable of all the selections stored in the shared key store.
        """
        for key, _ in self._state_dictionary.copy_content().items():  # type: ignore
            if key.startswith("selection."):
                yield self.get_selection(key)

    def get_selection(self, selection_id: str) -> RenderingSelection:
        """
        Get the selection with the given selection id, throwing a KeyError if
        it is not present. For the root selection, use the root_selection
        property.

        :param selection_id: The id of the selection
        :return: The selection if it is present
        """
        value = self._state_dictionary.copy_content()[selection_id]  # type: ignore
        return self._create_selection_from_dict(value)

    def _create_selection_from_dict(self, value) -> RenderingSelection:
        selection = RenderingSelection.from_dictionary(value)
        selection.updated.add_callback(self.update_selection)
        selection.removed.add_callback(self.remove_selection)
        return selection

    def _create_selection_from_id_and_name(
        self, selection_id: str, name: str
    ) -> RenderingSelection:
        selection = RenderingSelection(selection_id, name)
        selection.updated.add_callback(self.update_selection)
        selection.removed.add_callback(self.remove_selection)
        return selection
