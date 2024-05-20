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
        prev_yield = time.perf_counter()
        yield 0
        while True:
            time_since_yield = time.perf_counter() - prev_yield
            wait_duration = max(0.0, self.interval - time_since_yield)
            time.sleep(wait_duration)
            next_yield = time.perf_counter()
            yield next_yield - prev_yield
            prev_yield = next_yield


def yield_interval(interval: float):
    """
    Yield immediately and then every interval seconds, yielding the time in seconds that passed between yields.

    :param interval: Number of seconds to ensure between yields
    :yield: Number of seconds since last yielding
    """
    return VariableIntervalGenerator(interval).yield_interval()
