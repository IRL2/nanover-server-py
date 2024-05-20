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
        yield 0
        while True:
            wait_mixed(target - time.perf_counter())
            target += self.interval
            next_yield = time.perf_counter()
            yield next_yield - prev_yield
            prev_yield = next_yield


def yield_interval(interval: float):
    """
    Yield immediately and then roughly every interval seconds, yielding the time in seconds that passed between yields.

    :param interval: Number of seconds to ensure between yields
    :yield: Number of seconds since last yielding
    """
    return VariableIntervalGenerator(interval).yield_interval()


def wait_sleep(seconds: float):
    """
    Do nothing for a period of time using time.sleep. Imprecise but doesn't use much CPU.
    """
    target = time.perf_counter() + seconds
    duration = max(target - time.perf_counter(), 0)
    if duration > 0:
        time.sleep(duration)


def wait_busy(seconds: float):
    """
    Do nothing for a period of time by tightly looping. Precise but uses much CPU.
    """
    target = time.perf_counter() + seconds
    while time.perf_counter() < target:
        pass


def wait_mixed(seconds: float, sleep_error=.01):
    """
    Do nothing for a period of time by sleeping until close to the target then using CPU time.
    """
    target = time.perf_counter() + seconds
    while True:
        now = time.perf_counter()
        duration = max(target - now, 0)
        if duration == 0:
            return
        elif duration > sleep_error:
            time.sleep(duration * .5)
        else:
            while time.perf_counter() < target:
                pass

