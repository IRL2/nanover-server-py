import time


def yield_interval(interval):
    """
    Yield at a set interval, accounting for the time spent outside of this
    function.
    :param interval: Number of seconds to put between yields
    """
    last_yield = time.monotonic() - interval
    while True:
        time_since_yield = time.monotonic() - last_yield
        wait_duration = max(0, interval - time_since_yield)
        time.sleep(wait_duration)
        yield time.monotonic() - last_yield
        last_yield = time.monotonic()
