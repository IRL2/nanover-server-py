from contextlib import contextmanager

import pytest
from ase import Atoms
from nanover.imd import ImdServer, ImdClient, ImdStateWrapper
from nanover.state.state_dictionary import StateDictionary


@pytest.fixture
def imd_server():
    with ImdServer(address="localhost", port=0) as server:
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
def client_interaction(client: ImdClient, interaction):
    interaction_id = client.start_interaction()
    client.update_interaction(interaction_id, interaction)
    yield interaction_id
    client.stop_interaction(interaction_id)
