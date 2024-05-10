import time

import pytest
import itertools

from numpy import average

from nanover.utilities.timing import yield_interval

TIMING_TOLERANCE = 0.005  # 5ms
COMMON_INTERVALS = (1 / 10, 1 / 30, 1 / 60)


@pytest.mark.parametrize("interval", COMMON_INTERVALS)
def test_yield_interval(interval):
    """
    Test that yield_interval yields on average at the correct interval.
    """
    count = round(1 / interval)
    times = [
        time.monotonic() for _ in itertools.islice(yield_interval(interval), count)
    ]
    intervals = [times[i] - times[i - 1] for i in range(1, len(times))]
    assert average(intervals) == pytest.approx(interval, abs=TIMING_TOLERANCE)


@pytest.mark.parametrize("interval", COMMON_INTERVALS)
@pytest.mark.parametrize("delay", (0.5, 0.1, 0.01))
def test_yield_interval_delays(interval, delay):
    """
    Test that yield_interval does not yield any faster if more time is spent outside of it.
    """

    def do_delay():
        time.sleep(delay)
        return time.monotonic()

    # try not to test for longer than 1s but do at least 3 yields
    count = max(3, round(1 / (interval + delay)))
    times = [do_delay() for _ in itertools.islice(yield_interval(interval), count)]
    intervals = [times[i] - times[i - 1] for i in range(1, len(times))]
    assert average(intervals) == pytest.approx(interval + delay, abs=TIMING_TOLERANCE)


@pytest.mark.parametrize("interval", COMMON_INTERVALS)
@pytest.mark.parametrize("delay", (0.5, 0.1, 0.01))
def test_yield_interval_delays_dt(interval, delay):
    """
    Test that yield_interval yields only the time spent inside itself.
    """

    # try not to test for longer than 1s but do at least 3 yields
    count = max(3, round(1 / (interval + delay)))

    yield_times = []
    enter_times = []

    def monitor_delay(dt):
        yield_times.append(time.monotonic())
        time.sleep(delay)
        enter_times.append(time.monotonic())
        return dt

    enter_times.append(time.monotonic())
    reported_deltas = [
        monitor_delay(dt) for dt in itertools.islice(yield_interval(interval), count)
    ]
    measured_deltas = [yield_times[i] - enter_times[i] for i in range(count)]

    assert reported_deltas[1:] == pytest.approx(
        measured_deltas[1:], abs=TIMING_TOLERANCE
    )
