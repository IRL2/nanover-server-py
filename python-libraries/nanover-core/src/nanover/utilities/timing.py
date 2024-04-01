"""
Module providing time-dependent utility methods.
"""

import time
import os
from threading import RLock


# In NanoVer, the `time.sleep` function is consistently invoked between
# successive transmissions in a gRPC data stream to regulate the data transfer
# frequency. However, on Windows systems, the `time.sleep` function commonly
# fails to provide reliable delays, particularly for shorter intervals. To
# address this issue, a more dependable, low-level sleep function is employed.
# This solution will become less critical in the future once an application
# scheduler has been implemented.
if os.name == 'nt':
    import ctypes
    from ctypes.wintypes import LARGE_INTEGER

    def better_sleep(seconds):
        """Delay execution for a given number of seconds.

        :params seconds: number of seconds to sleep for

        """
        # Creates a waitable timer object and returns a handle to the object.
        handle = ctypes.windll.kernel32.CreateWaitableTimerExW(
            # LPSECURITY_ATTRIBUTES = null
            None,
            # LPCWSTR = null
            None,
            # CREATE_WAITABLE_TIMER_HIGH_RESOLUTION
            0x00000002,
            # EVENT_ALL_ACCESS
            0x1F0003
        )

        # Start the wait timer
        ctypes.windll.kernel32.SetWaitableTimer(
            # Waitable timer handle
            handle,
            # Time specified in 100 nanosecond intervals, must be negative to
            # indicate relative wait.
            ctypes.byref(LARGE_INTEGER(int(seconds * -10000000))),
            # Do not repeat timer
            0,
            # No callback routine required
            None,
            # No arguments required for the callback routine
            None,
            # Don't attempt to wake the computer from a suspended power state
            0,
        )

        # Wait for the timer to elapse
        ctypes.windll.kernel32.WaitForSingleObject(
            # Pointer to wait event
            handle,
            # Timeout for wait event (infinite)
            0xFFFFFFFF)

        # Clean up the wait timer
        ctypes.windll.kernel32.CancelWaitableTimer(handle)

else:
    # If not on Windows, then default to the standard time.sleep function
    better_sleep = time.sleep


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
            better_sleep(wait_duration)
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
        better_sleep(wait_duration)
        yield time.monotonic() - last_yield
        last_yield = time.monotonic()
