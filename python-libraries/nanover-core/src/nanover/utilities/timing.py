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
        yield 0
        while True:
            last_yield = time.monotonic()
            time.sleep(self.interval)
            yield time.monotonic() - last_yield


def yield_interval(interval: float):
    """
    Yield immediately then sleep for a specified time between each subsequent yield.

    :param interval: Number of seconds to sleep
    :yield: Number of seconds actually spent sleeping
    """
    yield 0
    while True:
        last_yield = time.monotonic()
        time.sleep(interval)
        yield time.monotonic() - last_yield
