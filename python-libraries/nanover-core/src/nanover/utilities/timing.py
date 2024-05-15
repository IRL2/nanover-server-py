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
        target = time.perf_counter()
        prev_yield = target
        while True:
            while time.perf_counter() < target:
                time.sleep(max(0.0, (target - time.perf_counter()) * 0.5))
            target += self.interval
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


def sleep_precise(duration: float):
    prev_time = time.perf_counter()
    next_time = prev_time + duration
    while time.perf_counter() < next_time:
        sleep_time = next_time - time.perf_counter()
        time.sleep(sleep_time * 0.5)
