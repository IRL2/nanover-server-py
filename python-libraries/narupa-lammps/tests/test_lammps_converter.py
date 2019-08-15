import pytest
from narupa.lammps import LammpsHook
from narupa.lammps import DummyLammps
from narupa.trajectory.frame_data import POSITIONS, ELEMENTS
from narupa.trajectory import FrameServer, FrameData

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
    assert len(frame_data.raw.arrays[POSITIONS].float_values.values) == 9
    del h
    #FrameServer.close()

def test_elements_lammps_atoms():
    """
    Checks that the dimensionality of the position array is correctly returned for dummy LAMMPS
    """
    h = LammpsHook()
    dummy = DummyLammps(3)
    dummy.default_atoms = 3
    frame_data = FrameData()
    atom_type, masses = h.gather_lammps_particle_types(dummy)
    frame_data.arrays[ELEMENTS] = atom_type
    assert frame_data.raw.arrays[ELEMENTS].index_values.values == [1, 1, 1]
