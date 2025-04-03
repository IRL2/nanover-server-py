import time
from signal import signal, SIGINT
from contextlib import contextmanager


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
        self._cancelled = True

    def wait_cancellation(self, interval=0.1):
        """
        Sleep until this token is cancelled.

        :param interval: the interval in seconds between waking up to check cancellation.
        """
        while not self._cancelled:
            time.sleep(interval)
