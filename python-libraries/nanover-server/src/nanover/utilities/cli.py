import time
from signal import signal, SIGINT
from contextlib import contextmanager

from nanover.utilities.event import Event


@contextmanager
def suppress_keyboard_interrupt_as_cancellation():
    """
    Context manager that suppresses KeyboardInterrupt and instead yields a cancellation token that can be polled to
    check if a keyboard interrupt has occurred during the lifetime of the context.
    """
    token = CancellationToken()

    prev_handler = signal(SIGINT, lambda _, __: token.cancel())
    try:
        yield token
    finally:
        signal(SIGINT, prev_handler)


class CancellationToken:
    _cancelled = False
    _on_cancellation = Event()

    def subscribe_cancellation(self, callback):
        if self._cancelled:
            callback()
        else:
            self._on_cancellation.add_callback(callback)

    @property
    def is_cancelled(self):
        """
        Has this token been cancelled?
        """
        return self._cancelled

    def cancel(self):
        """
        Cancel this token.
        """
        if not self._cancelled:
            self._cancelled = True
            self._on_cancellation.invoke()

    def wait_cancellation(self, interval=0.1):
        """
        Sleep until this token is cancelled.

        :param interval: the interval in seconds between waking up to check cancellation.
        """
        while not self._cancelled:
            time.sleep(interval)
