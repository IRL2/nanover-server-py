import time
from typing import Callable


def assert_equal_soon(a: Callable, b: Callable, interval=0.1, timeout=1.0):
    __tracebackhide__ = True
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline and a() != b():
        time.sleep(interval)
    assert a() == b()


def assert_in_soon(a: Callable, b: Callable, interval=0.1, timeout=1.0):
    __tracebackhide__ = True
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline and a() not in b():
        time.sleep(interval)
    assert a() in b()


def assert_not_in_soon(a: Callable, b: Callable, interval=0.1, timeout=1.0):
    __tracebackhide__ = True
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline and a() in b():
        time.sleep(interval)
    assert a() not in b()
