from nanover.utilities.event import Event
from unittest.mock import Mock


def test_invoke_nocallbacks():
    event = Event()
    event.invoke("argument")


def test_invoke_callback():
    event = Event()

    def callback():
        callback.called = True

    callback.called = False

    event.add_callback(callback)

    assert callback.called is False

    event.invoke()

    assert callback.called is True


def test_invoke_callback_then_remove():
    event = Event()

    def callback():
        callback.called += 1

    callback.called = 0

    event.add_callback(callback)

    assert callback.called == 0

    event.invoke()

    assert callback.called == 1

    event.remove_callback(callback)

    event.invoke()

    assert callback.called == 1


def test_exceptions_suppressed():
    """Test exceptions in callbacks are suppressed and that all callbacks are called regardless."""
    event = Event()
    callback1 = Mock(side_effect=Exception)
    callback2 = Mock()
    event.add_callback(callback1)
    event.add_callback(callback2)
    event.invoke()
    assert callback1.called and callback2.called
