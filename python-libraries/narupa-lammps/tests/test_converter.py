import pytest
import narupa.lammps.hook as nlh
from narupa.lammps import LammpsHook

@pytest.fixture
def simple_atom_lammps_frame():
    n_atoms = 3
    data_array = nlh.manipulate_dummy_array("x", n_atoms)
    return data_array


def test_topology_lammps_atoms(simple_atom_lammps_frame):
    h = LammpsHook()
    frame_data = h.lammps_to_frame_data(simple_atom_lammps_frame, positions=True, topology=False)
    assert len(frame_data.arrays['atom.position'].float_values.values) == 9


