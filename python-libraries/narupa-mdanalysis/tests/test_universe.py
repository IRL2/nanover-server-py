import os
import pytest
import numpy as np
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
USER_FORCES_TRAJ = os.path.join(
    os.path.dirname(os.path.realpath(__file__)),
    "hello_force.traj",
)
REFERENCE_PDB = os.path.join(
    os.path.dirname(os.path.realpath(__file__)),
    "17-ala.pdb",
)


@pytest.fixture
def single_topology_universe():
    return mda.Universe(
        SINGLE_TOPOLOGY_TRAJ,
        format=NarupaReader,
        topology_format=NarupaParser,
    )


@pytest.fixture
def multi_topology_universe():
    return mda.Universe(
        MULTI_TOPOLOGY_TRAJ,
        format=NarupaReader,
        topology_format=NarupaParser,
    )


@pytest.fixture
def user_forces_universe():
    return mda.Universe(
        USER_FORCES_TRAJ,
        format=NarupaReader,
        topology_format=NarupaParser,
    )


@pytest.fixture
def reference_topology_universe():
    return mda.Universe(REFERENCE_PDB)


def test_n_frames(single_topology_universe):
    assert single_topology_universe.trajectory.n_frames == 623


def test_n_frames_multi_topology(multi_topology_universe):
    assert multi_topology_universe.trajectory.n_frames == 2024


def test_n_atoms(single_topology_universe):
    assert len(single_topology_universe.atoms) == 173


@pytest.mark.parametrize(
    "attribute",
    (
        "names",
        "resnames",
        "resindices",
        "chainIDs",
        "segids",
        "elements",
        "types",
        "segindices",
    ),
)
def test_topology_attributes(
    attribute, single_topology_universe, reference_topology_universe
):
    actual_attribute = getattr(single_topology_universe.atoms, attribute)
    expected_attribute = getattr(reference_topology_universe.atoms, attribute)
    assert all(actual_attribute == expected_attribute)


def test_user_forces_all_frames(user_forces_universe):
    """
    When forces are available, they are for all the frames.
    """
    assert all(
        "user_forces" in ts.data
        for ts in user_forces_universe.trajectory
    )


def test_user_forces_shape(user_forces_universe):
    """
    When forces are available, the shape of the array is correct.
    """
    shape = (len(user_forces_universe.atoms), 3)
    assert all(
        ts.data["user_forces"].shape == shape
        for ts in user_forces_universe.trajectory
    )


def test_user_forces(user_forces_universe):
    """
    There are the expected number of frames with non zero forces.

    This is a regression test, the expected number is obtained from running the
    code.
    """
    expected = 414
    actual = sum(
        np.any(ts.data['user_forces'])
        for ts in user_forces_universe.trajectory
    )
    assert actual == expected
