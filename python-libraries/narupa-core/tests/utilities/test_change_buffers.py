# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.

import pytest

from narupa.utilities.change_buffers import (
    DictionaryChangeBuffer, DictionaryChangeMultiView, ObjectFrozenError,
)


@pytest.fixture
def change_buffer():
    return DictionaryChangeBuffer()


@pytest.fixture
def change_multiview():
    return DictionaryChangeMultiView()


def test_buffer_flush_reflects_changes(change_buffer):
    """
    Test that flushing reflects the previous update.
    """
    change_buffer.update({"hello": "test"})
    changes, removals = change_buffer.flush_changed_blocking()
    assert changes["hello"] == "test"


def test_buffer_flush_reflects_removal(change_buffer):
    """
    Test that flushing reflects the previous key removal.
    """
    change_buffer.update({"hello": "test"})
    change_buffer.update(removals=["hello"])
    changes, removals = change_buffer.flush_changed_blocking()
    assert removals == set(["hello"])


def test_buffer_flush_empties_changes(change_buffer):
    """
    Test that flushing empties the buffer of changes.
    """
    change_buffer.update({"hello": "test"})
    assert change_buffer._pending_changes
    change_buffer.flush_changed_blocking()
    assert not change_buffer._pending_changes


def test_buffer_flush_empties_removals(change_buffer):
    """
    Test that flushing empties the buffer of removals.
    """
    change_buffer.update({"hello": "test"})
    change_buffer.update(removals=["hello"])
    assert change_buffer._pending_removals
    change_buffer.flush_changed_blocking()
    assert not change_buffer._pending_removals


def test_buffer_flush_merges_updates(change_buffer):
    """
    Test that flushing after two updates gives a single combined update.
    """
    change_buffer.update({"hello": "test"})
    change_buffer.update({"foo": "bar"})
    changes, removals = change_buffer.flush_changed_blocking()
    assert changes["hello"] == "test" and changes["foo"] == "bar"


def test_buffer_flush_merges_removals(change_buffer):
    """
    Test that flushing after two removals gives a single combined removal.
    """
    change_buffer.update({"hello": "test", "foo": "bar"})
    change_buffer.update(removals=["hello"])
    change_buffer.update(removals=["foo"])
    changes, removals = change_buffer.flush_changed_blocking()
    assert removals == set(["hello", "foo"])


def test_buffer_flush_merges_same_change_key(change_buffer):
    """
    Test that flushing after two updates of the same key gives a single latest
    value.
    """
    change_buffer.update({"hello": "test"})
    change_buffer.update({"hello": "bar"})
    changes, removals = change_buffer.flush_changed_blocking()
    assert changes["hello"] == "bar"


def test_change_then_removal_discards_change(change_buffer):
    """
    Test that flushing after adding a key then removing the key gives no change.
    """
    change_buffer.update({"hello": "test"})
    change_buffer.update(removals=["hello"])
    changes, removals = change_buffer.flush_changed_blocking()
    assert not changes


def test_frozen_buffer_cant_update(change_buffer):
    """
    Test that attempting to update after freezing the buffer raises the correct
    exception.
    """
    change_buffer.freeze()
    with pytest.raises(ObjectFrozenError):
        change_buffer.update({"hello": "test"})


@pytest.mark.timeout(1)
def test_frozen_empty_buffer_cant_flush(change_buffer):
    """
    Test that flushing an empty buffer after freezing it raises the correct
    exception.
    """
    change_buffer.freeze()
    with pytest.raises(ObjectFrozenError):
        change_buffer.flush_changed_blocking()


def test_frozen_buffer_update_ignored(change_buffer):
    """
    Test that a failed update after freezing a buffer does not affect the
    unflushed changes.
    """
    change_buffer.update({"hello": "test"})
    change_buffer.freeze()
    try:
        change_buffer.update({"foo": "bar"})
    except ObjectFrozenError:
        pass
    changes, removals = change_buffer.flush_changed_blocking()
    assert changes["hello"] == "test" and "foo" not in changes


def test_frozen_multiview_cant_update(change_multiview):
    """
    Test that attempting to update a frozen multiview raises the correct
    exception.
    """
    change_multiview.freeze()
    with pytest.raises(ObjectFrozenError):
        change_multiview.update({"hello": "test"})


@pytest.mark.timeout(1)
def test_frozen_multiview_view_gives_last_values_and_no_removals(change_multiview):
    """
    Test that views can still be created on a frozen multiview but that they
    only provide the initial values and then raise the correct exception on
    subsequent flushes.
    """
    change_multiview.update({"hello": "test", "foo": "bar"})
    change_multiview.update(removals={"foo"})
    change_multiview.freeze()
    with change_multiview.create_view() as view:
        changes, removals = view.flush_changed_blocking()
        assert changes["hello"] == "test"
        assert "foo" not in changes
        assert not removals
        with pytest.raises(ObjectFrozenError):
            view.flush_changed_blocking()


@pytest.mark.timeout(1)
def test_frozen_multiview_subscribe_gives_last_values_and_no_removals(change_multiview):
    """
    Test that subscribing a frozen multiview provides the initial values and
    then ends.
    """
    change_multiview.update({"hello": "test"})
    change_multiview.freeze()
    for changes, removals in change_multiview.subscribe_changes():
        assert changes["hello"] == "test" and not removals
