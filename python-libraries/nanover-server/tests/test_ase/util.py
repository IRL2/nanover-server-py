from contextlib import contextmanager

import pytest
from ase import Atoms

from nanover.app import NanoverImdApplication
from nanover.essd.utils import get_broadcastable_ip
from nanover.imd import ImdStateWrapper
from nanover.state.state_dictionary import StateDictionary


@pytest.fixture
def app_server():
    with NanoverImdApplication.basic_server(
        address=get_broadcastable_ip(), port=0
    ) as server:
        yield server


@pytest.fixture
def state_wrapper():
    return ImdStateWrapper(StateDictionary())


def co_atoms():
    d = 1.1
    co = Atoms("CO", positions=[(0, 0, 0), (0, 0, d)], cell=[2, 2, 2], pbc=[1, 1, 1])
    return co


def c_atoms():
    c = Atoms("C", positions=[(0, 0, 0)], cell=[2, 2, 2], pbc=[1, 1, 1])
    return c


@contextmanager
def client_interaction(client, interaction):
    interaction_id = client.start_interaction()
    client.update_interaction(interaction_id, interaction)
    yield interaction_id
    client.stop_interaction(interaction_id)
