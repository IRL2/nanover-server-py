import pytest
from narupa.core.key_lockable_map import ResourceLockedException
from narupa.multiplayer.change_buffers import DictionaryChange
from narupa.state import StateDictionary


ACCESS_TOKEN_1 = object()
ACCESS_TOKEN_2 = object()

INITIAL_STATE = {
    'hello': 100,
    'test': {'baby': 'yoda'},
}


@pytest.fixture
def state_dictionary():
    state_dictionary = StateDictionary()
    change = DictionaryChange(INITIAL_STATE, set())
    state_dictionary.update_state(None, change)
    return state_dictionary


def test_initial_state(state_dictionary):
    """
    Test that the initial state of the dictionary gets set.
    """
    assert state_dictionary.content == INITIAL_STATE


def test_update_unlocked(state_dictionary):
    """
    Test that unlocked keys can be changed and removed.
    """
    update = DictionaryChange({'hello': 50}, {'test'})
    state_dictionary.update_state(ACCESS_TOKEN_1, update)

    assert state_dictionary.content == {'hello': 50, }


def test_partial_lock_atomic(state_dictionary):
    """
    Test that an update attempt has no effect if the whole update cannot be
    made.
    """
    state_dictionary.update_locks(ACCESS_TOKEN_2, {'hello': 10}, set())
    update = DictionaryChange({'hello': 50, 'goodbye': 50}, set())

    with pytest.raises(ResourceLockedException):
        state_dictionary.update_state(ACCESS_TOKEN_1, update)

    assert state_dictionary.content == INITIAL_STATE


def test_unheld_releases_ignored(state_dictionary):
    """
    Test that attempting to release locks on keys which are not locked with this
    access token will not prevent the update from occurring.
    """
    update = DictionaryChange({'hello': 50}, {'goodbye'})
    state_dictionary.update_state(ACCESS_TOKEN_1, update)
    assert state_dictionary.content['hello'] == 50
