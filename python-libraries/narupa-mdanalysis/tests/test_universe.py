import os
import pytest
import MDAnalysis as mda
from narupa.mdanalysis import NarupaParser, NarupaReader

SINGLE_TOPOLOGY_TRAJ = os.path.join(
    os.path.dirname(os.path.realpath(__file__)),
    "hello.traj",
)
MULTI_TOPOLOGY_TRAJ = os.path.join(
    os.path.dirname(os.path.realpath(__file__)),
    "hello_multi.traj",
)
REFERENCE_PDB = os.path.join(
    os.path.dirname(os.path.realpath(__file__)),
    "17-ala.pdb",
)


@pytest.fixture
def single_topology_universe():
    return mda.Universe(
        SINGLE_TOPOLOGY_TRAJ,
        SINGLE_TOPOLOGY_TRAJ,
        format=NarupaReader,
        topology_format=NarupaParser,
    )


@pytest.fixture
def multi_topology_universe():
    return mda.Universe(
        MULTI_TOPOLOGY_TRAJ,
        MULTI_TOPOLOGY_TRAJ,
        format=NarupaReader,
        topology_format=NarupaParser,
    )


@pytest.fixture
def reference_topology_universe():
    return mda.Universe(REFERENCE_PDB)


def test_n_frames(single_topology_universe):
    assert single_topology_universe.trajectory.n_frames == 622


def test_n_frames_multi_topology(multi_topology_universe):
    assert multi_topology_universe.trajectory.n_frames == 2023


def test_n_atoms(single_topology_universe):
    assert len(single_topology_universe.atoms) == 173


@pytest.mark.parametrize('attribute', (
    'names', 'resnames', 'resindices', 'chainIDs', 'segids', 'elements', 'types', 'segindices', 
))
def test_topology_attributes(attribute, single_topology_universe, reference_topology_universe):
    actual_attribute = getattr(single_topology_universe.atoms, attribute)
    expected_attribute = getattr(reference_topology_universe.atoms, attribute)
    assert all(actual_attribute == expected_attribute)

