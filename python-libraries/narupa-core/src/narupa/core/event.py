class Event:
    """
    A class which can have callback added and removed using += and -=, and invokes them when called.
    """

    def __init__(self):
        self._callbacks = []

    def __iadd__(self, callback):
        self._callbacks.append(callback)
        return self

    def __isub__(self, callback):
        self._callbacks.remove(callback)
        return self

    def __call__(self, *args, **kwargs):
        for callback in self._callbacks:
            callback(*args, **kwargs)