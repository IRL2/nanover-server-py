import pytest
import numpy as np
from narupa.lammps import LammpsHook
from narupa.lammps import DummyLammps
from narupa.trajectory.frame_data import POSITIONS, ELEMENTS
from narupa.trajectory import FrameData


@pytest.fixture
def simple_atom_lammps_frame():
    dummy = DummyLammps()
    data_array = dummy.gather_atoms("x", "", "")
    return data_array


def test_length_lammps_atoms(simple_atom_lammps_frame):
    """
    Checks that the dimensionality of the position array is correctly returned for dummy LAMMPS
    """
    h = LammpsHook()
    h.default_atoms = 3

    frame_data = FrameData()
    h.distance_factor = 1.0
    h.lammps_positions_to_frame_data(frame_data, simple_atom_lammps_frame)
    h.close()
    assert len(frame_data.raw.arrays[POSITIONS].float_values.values) == 9


def test_elements_lammps_atoms():
    """
    Checks that the dimensionality of the position array is correctly returned for dummy LAMMPS
    """
    # Instantiate classes
    h = LammpsHook()
    dummy = DummyLammps(3)
    # Check that get_atoms works
    dummy.default_atoms = 3
    assert dummy.get_natoms() == 3

    frame_data = FrameData()
    atom_type, masses = h.gather_lammps_particle_types(dummy)
    frame_data.arrays[ELEMENTS] = atom_type
    h.close()
    assert frame_data.raw.arrays[ELEMENTS].index_values.values == [1, 1, 1]


def test_forces_lammps_atoms():
    """
    Checks that the c_type conversion of forces for dummy LAMMPS

    """
    h = LammpsHook()
    h.force_factor = 1.0
    h.n_atoms = 3
    dummy = DummyLammps()
    dummy.default_atoms = 3
    # Collect empty  lammps c_type array
    data_array = h.manipulate_lammps_array("f", dummy)

    # Set external forces as numpy array
    temp_forces = np.array([1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0])
    h.add_interaction_to_ctype(temp_forces, data_array)
    # Convert back to numpy for assert
    final_forces = np.ctypeslib.as_array(data_array)
    h.close()
    assert final_forces.all() == np.array([1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0]).all()


def test_main_hook():
    """
    Checks the main hook routine work with all the above commands
    """
    h = LammpsHook()
    h.default_atoms = 3
    h.lammps_hook()
    # Get first three positions and convert to list floats
    positions = h.frame_data.raw.arrays[POSITIONS].float_values.values[0:3]
    positions = [float("%.1f" % x) for x in positions]
    h.close()
    assert positions == [0.0, 0.1, 0.2]
