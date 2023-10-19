import pytest
import numpy as np
from narupa.lammps import LammpsImd
from narupa.lammps.mock import MockLammps
from narupa.trajectory.frame_data import PARTICLE_POSITIONS, PARTICLE_ELEMENTS
from narupa.trajectory import FrameData
import sys

pytestmark = pytest.mark.skipif(
    sys.platform == "win32",
    reason=(
        "These tests can break the windows test runner on github. "
        "This is tracked in issue #33: "
        "<https://github.com/IRL2/narupa-protocol/issues/33>."
    ),
)


@pytest.fixture
def simple_atom_lammps_frame():
    mock = MockLammps()
    data_array = mock.gather_atoms("x", "", "")
    return data_array


@pytest.fixture
def lammps_hook():
    hook = LammpsImd()
    yield hook
    hook.close()


def test_length_lammps_atoms(simple_atom_lammps_frame, lammps_hook):
    """
    Checks that the dimensionality of the position array is correctly returned for mock LAMMPS
    """

    lammps_hook.n_atoms_in_dummy = 3

    frame_data = FrameData()
    lammps_hook.distance_factor = 1.0
    lammps_hook._lammps_positions_to_frame_data(frame_data, simple_atom_lammps_frame)
    assert len(frame_data.raw.arrays[PARTICLE_POSITIONS].float_values.values) == 9


def test_get_atoms(lammps_hook):
    """
    Checks that the atoms are correctly set withing Dummy_Lammps
    """
    # Instantiate classes
    mock = MockLammps(3)
    # Check that get_atoms works
    mock.n_atoms_in_dummy = 3
    assert mock.get_natoms() == 3


def test_elements_lammps_atoms(lammps_hook):
    """
    Checks that the dimensionality of the position array is correctly returned for mock LAMMPS
    """
    # Instantiate classes
    mock = MockLammps(3)
    # Check that get_atoms works
    mock.n_atoms_in_dummy = 3
    frame_data = FrameData()
    atom_type, masses = lammps_hook._gather_lammps_particle_types(mock)
    frame_data.arrays[PARTICLE_ELEMENTS] = atom_type
    assert frame_data.raw.arrays[PARTICLE_ELEMENTS].index_values.values == [1, 1, 1]


def test_forces_lammps_atoms(lammps_hook):
    """
    Checks that the c_type conversion of forces for mock LAMMPS

    """
    lammps_hook.force_factor = 1.0
    lammps_hook.n_atoms = 3
    mock = MockLammps()
    mock.n_atomms_in_dummy = 3
    # Collect empty  lammps c_type array
    data_array = lammps_hook._gather_lammps_array("f", mock)

    # Set external forces as numpy array
    temp_forces = np.array([1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0])
    lammps_hook._add_interaction_to_ctype(temp_forces, data_array)
    # Convert back to numpy for assert
    final_forces = np.ctypeslib.as_array(data_array)
    assert (
        final_forces.all()
        == np.array([1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0]).all()
    )


def test_main_hook(lammps_hook):
    """
    Checks the main hook routine work with all the above commands
    """
    lammps_hook.default_atoms = 3
    lammps_hook.lammps_hook()
    # Get first three positions and convert to list floats
    positions = lammps_hook.frame_data.raw.arrays[
        PARTICLE_POSITIONS
    ].float_values.values[0:3]
    positions = [float("%.1f" % x) for x in positions]
    assert positions == [0.0, 0.1, 0.2]
