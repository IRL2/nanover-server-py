import time

import numpy
import pytest
import itertools

from numpy import average

from nanover.utilities.timing import yield_interval

TIMING_TOLERANCE = 0.005  # 5ms
COMMON_INTERVALS = (1 / 10, 1 / 30, 1 / 60)


@pytest.mark.serial  # we want accurate timing so run without any parallel load
@pytest.mark.parametrize("interval", COMMON_INTERVALS)
@pytest.mark.parametrize("work_factor", (0.75, 0.5, 0.25, 0))
def test_yield_interval(interval, work_factor):
    """
    Test that yield_interval yields on average at the correct interval when time is spent between resuming iteration.
    """

    times = []
    count = round(1 / interval)

    for dt in itertools.islice(yield_interval(interval), count):
        times.append(time.perf_counter())
        time.sleep(interval * work_factor)

    intervals = numpy.diff(times)
    assert average(intervals) == pytest.approx(interval, abs=TIMING_TOLERANCE)


@pytest.mark.parametrize("interval", COMMON_INTERVALS)
@pytest.mark.parametrize("work_factor", (0.75, 0.5, 0.25, 0))
def test_yield_interval_dt(interval, work_factor):
    """
    Test that yield_interval yields the actual time between yields.
    """

    times = []
    reported_deltas = []
    count = round(1 / interval)

    times.append(time.perf_counter())
    for dt in itertools.islice(yield_interval(interval), count):
        reported_deltas.append(dt)
        times.append(time.perf_counter())
        time.sleep(interval * work_factor)

    measured_deltas = numpy.diff(times)
    assert reported_deltas == pytest.approx(measured_deltas, abs=TIMING_TOLERANCE)
