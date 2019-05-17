import pytest
from narupa.lammps import LammpsHook
from narupa.lammps import DummyLammps
from narupa.trajectory.frame_data import POSITIONS
from narupa.protocol.trajectory import FrameData
from narupa.trajectory import FrameServer, FrameData

@pytest.fixture
def simple_atom_lammps_frame():
    n_atoms = 3
    dummy = DummyLammps()
    data_array = dummy.manipulate_dummy_array(n_atoms)
    return data_array


def test_topology_lammps_atoms(simple_atom_lammps_frame):
    h = LammpsHook()
    frame_data_obj = FrameData()
    frame_data = h.lammps_positions_to_frame_data(frame_data_obj, simple_atom_lammps_frame)
    assert len(frame_data.raw.arrays[POSITIONS].float_values.values) == 9


