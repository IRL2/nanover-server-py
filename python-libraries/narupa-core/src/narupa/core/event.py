from typing import Callable, TypeVar

class Event:
    """
    A class which can have callback added and removed using += and -=, and invokes them when called.
    """

    def __init__(self):
        self._callbacks = []

    def add_callback(self, callback: Callable[..., None]):
        """
        Add a callback to this event, which will be invoked everytime this event is invoked.

        :param callback: The callback to be called when this event is triggered
        """
        self._callbacks.append(callback)

    def remove_callback(self, callback: Callable[..., None]):
        """
        Remove a callback from this event.

        :param callback: The callback to be removed from this event's callbacks
        """
        self._callbacks.remove(callback)

    def __call__(self, *args, **kwargs):
        for callback in self._callbacks:
            callback(*args, **kwargs)