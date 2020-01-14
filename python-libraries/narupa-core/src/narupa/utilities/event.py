from typing import Callable


class Event:
    """
    A class which stores a set of callbacks, which are invoked when an event is published.
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

    def invoke(self, *args, **kwargs):
        """
        Invoke the callbacks associated with this event with the provided arguments.

        :param args: Positional arguments for the event, passed on to each callback.
        :param kwargs: Keywords arguments for the event, passed on to each callback.
        """
        for callback in self._callbacks:
            callback(*args, **kwargs)