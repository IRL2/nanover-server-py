# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.

import pytest

from narupa.multiplayer.dictionary_change_buffer import DictionaryChangeBuffer, DictionaryChangeMultiView, ObjectClosedException


@pytest.fixture
def change_buffer():
    yield DictionaryChangeBuffer()


@pytest.fixture
def change_multiview():
    yield DictionaryChangeMultiView()


def test_buffer_flush_reflects_update(change_buffer):
    """Test that flushing reflects the previous update."""
    change_buffer.update({"hello": "test"})
    changes = change_buffer.flush_changed_blocking()
    assert changes["hello"] == "test"


def test_buffer_flush_merges_updates(change_buffer):
    """Test that flushing after two updates gives a single combined update."""
    change_buffer.update({"hello": "test"})
    change_buffer.update({"foo": "bar"})
    changes = change_buffer.flush_changed_blocking()
    assert changes["hello"] == "test" and changes["foo"] == "bar"


def test_closed_buffer_cant_update(change_buffer):
    """Test that attempting to update after closing the buffer raises the
    correct exception."""
    change_buffer.close()
    with pytest.raises(ObjectClosedException):
        change_buffer.update({"hello": "test"})


@pytest.mark.timeout(1)
def test_closed_empty_buffer_cant_flush(change_buffer):
    """Test that flushing an empty buffer after closing it raises the correct
    exception."""
    change_buffer.close()
    with pytest.raises(ObjectClosedException):
        change_buffer.flush_changed_blocking()


def test_closed_buffer_update_ignored(change_buffer):
    """Test that a failed update after closing a buffer does not affect the
    unflushed changes."""
    change_buffer.update({"hello": "test"})
    change_buffer.close()
    try:
        change_buffer.update({"foo": "bar"})
    except ObjectClosedException:
        pass
    changes = change_buffer.flush_changed_blocking()
    assert changes["hello"] == "test" and "foo" not in changes


def test_closed_multiview_cant_update(change_multiview):
    """Test that attempting to update a close multiview raises the correct
    exception."""
    change_multiview.close()
    with pytest.raises(ObjectClosedException):
        change_multiview.update({"hello": "test"})


@pytest.mark.timeout(1)
def test_closed_multiview_view_gives_last_values(change_multiview):
    """Test that views can still be created on a closed multiview but that
    they only provide the initial values and then raise the correct exception
    on subsequent flushes."""
    change_multiview.update({"hello": "test"})
    change_multiview.close()
    view = change_multiview.create_view()
    changes = view.flush_changed_blocking()
    assert changes["hello"] == "test"
    with pytest.raises(ObjectClosedException):
        view.flush_changed_blocking()


@pytest.mark.timeout(1)
def test_closed_multiview_subscribe_gives_last_values(change_multiview):
    """Test that subscribing a closed multiview provides the initial values and
    then ends."""
    change_multiview.update({"hello": "test"})
    change_multiview.close()
    for changes in change_multiview.subscribe_updates():
        assert changes["hello"] == "test"
