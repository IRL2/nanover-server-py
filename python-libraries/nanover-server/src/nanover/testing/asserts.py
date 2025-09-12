import time
from contextlib import suppress
from typing import Callable


DEFAULT_INTERVAL = 0.1
DEFAULT_TIMEOUT = 0.25


def assert_true_soon(
    p: Callable, *, interval=DEFAULT_INTERVAL, timeout=DEFAULT_TIMEOUT
):
    __tracebackhide__ = True  # hide this function in the test traceback
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        with suppress(Exception):
            if p():
                break
        time.sleep(interval)
    assert p()


def assert_equal_soon(
    a: Callable, b: Callable, *, interval=DEFAULT_INTERVAL, timeout=DEFAULT_TIMEOUT
):
    __tracebackhide__ = True  # hide this function in the test traceback
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        with suppress(Exception):
            if a() == b():
                break
        time.sleep(interval)
    assert a() == b()


def assert_in_soon(
    a: Callable, b: Callable, *, interval=DEFAULT_INTERVAL, timeout=DEFAULT_TIMEOUT
):
    __tracebackhide__ = True  # hide this function in the test traceback
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        with suppress(Exception):
            if a() in b():
                break
        time.sleep(interval)
    assert a() in b()


def assert_not_in_soon(
    a: Callable, b: Callable, *, interval=DEFAULT_INTERVAL, timeout=DEFAULT_TIMEOUT
):
    __tracebackhide__ = True  # hide this function in the test traceback
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        with suppress(Exception):
            if a() not in b():
                break
        time.sleep(interval)
    assert a() not in b()
