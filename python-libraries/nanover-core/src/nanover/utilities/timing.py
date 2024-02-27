"""
Module providing time-dependent utility methods.
"""
import time
from threading import RLock


class VariableIntervalGenerator:
    def __init__(self, default_interval):
        self._interval_lock = RLock()
        self.interval = default_interval

    @property
    def interval(self):
        with self._interval_lock:
            return self._interval

    @interval.setter
    def interval(self, value):
        with self._interval_lock:
            self._interval = value

    def yield_interval(self):
        last_yield = time.monotonic() - self.interval
        while True:
            time_since_yield = time.monotonic() - last_yield
            wait_duration = max(0.0, self.interval - time_since_yield)
            time.sleep(wait_duration)
            yield time.monotonic() - last_yield
            last_yield = time.monotonic()


def yield_interval(interval: float):
    """
    Yield at a set interval, accounting for the time spent outside of this
    function.

    :param interval: Number of seconds to put between yields
    """
    last_yield = time.monotonic() - interval
    while True:
        time_since_yield = time.monotonic() - last_yield
        wait_duration = max(0.0, interval - time_since_yield)
        time.sleep(wait_duration)
        yield time.monotonic() - last_yield
        last_yield = time.monotonic()
