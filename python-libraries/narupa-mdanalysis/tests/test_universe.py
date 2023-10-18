import os
import pytest
import itertools
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
FORCES_TRAJ = os.path.join(
    os.path.dirname(os.path.realpath(__file__)),
    "hello_all_forces.traj",
)
VELOCITIES_TRAJ = os.path.join(
    os.path.dirname(os.path.realpath(__file__)),
    "hello_vel.traj",
)
VELOCITIES_FORCES_TRAJ = os.path.join(
    os.path.dirname(os.path.realpath(__file__)),
    "hello_vel_force.traj",
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


@pytest.fixture(
    params=(
        (True, True),
        (False, False),
        (True, False),
        (False, True),
    )
)
def feature_universe_and_features(request):
    with_velocities, with_forces = request.param
    possible_paths = {
        (True, True): VELOCITIES_FORCES_TRAJ,
        (False, False): SINGLE_TOPOLOGY_TRAJ,
        (True, False): VELOCITIES_TRAJ,
        (False, True): FORCES_TRAJ,
    }
    file_path = possible_paths[(with_velocities, with_forces)]
    return (
        mda.Universe(
            file_path,
            format=NarupaReader,
            topology_format=NarupaParser,
        ),
        with_velocities,
        with_forces,
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
    assert all("user_forces" in ts.data for ts in user_forces_universe.trajectory)


def test_user_forces_shape(user_forces_universe):
    """
    When forces are available, the shape of the array is correct.
    """
    shape = (len(user_forces_universe.atoms), 3)
    assert all(
        ts.data["user_forces"].shape == shape for ts in user_forces_universe.trajectory
    )


def test_user_forces(user_forces_universe):
    """
    There are the expected number of frames with non zero forces.

    This is a regression test, the expected number is obtained from running the
    code.
    """
    expected = 414
    actual = sum(
        np.any(ts.data["user_forces"]) for ts in user_forces_universe.trajectory
    )
    assert actual == expected


def ts_has_velocities(ts):
    try:
        ts.velocities
    except mda.exceptions.NoDataError:
        return False
    else:
        return True


def ts_has_forces(ts):
    try:
        ts.forces
    except mda.exceptions.NoDataError:
        return False
    else:
        return True


def test_velocities(feature_universe_and_features):
    """
    Velocities are optional in the recording.
    """
    universe, with_velocities, with_forces = feature_universe_and_features
    all_frame_have_velocities = all(ts_has_velocities(ts) for ts in universe.trajectory)
    assert with_velocities == all_frame_have_velocities


def test_forces(feature_universe_and_features):
    """
    Forces are optional in the recording.
    """
    universe, with_velocities, with_forces = feature_universe_and_features
    all_frame_have_forces = all(ts_has_forces(ts) for ts in universe.trajectory)
    assert with_forces == all_frame_have_forces
