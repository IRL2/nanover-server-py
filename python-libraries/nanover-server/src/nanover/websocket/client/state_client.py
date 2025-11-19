from typing import Any, ContextManager

from nanover.core.app_server import StateService
from nanover.utilities.change_buffers import DictionaryChange
from nanover.utilities.state_dictionary import StateDictionary
from nanover.websocket.client.base_client import WebsocketClient


class StateClient(WebsocketClient, StateService):
    """
    Mixin of methods for implementing StateService
    """

    def lock_state(self) -> ContextManager[dict[str, Any]]:
        return self._state_dictionary.lock_content()

    def copy_state(self) -> dict[str, Any]:
        return self.state_dictionary.copy_content()

    def update_state(self, change: DictionaryChange) -> None:
        self.state_dictionary.update_state(None, change)

    def clear_locks(self) -> None:
        raise NotImplementedError()

    @property
    def state_dictionary(self) -> StateDictionary:
        return self._state_dictionary
