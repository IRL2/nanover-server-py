import pytest
from narupa.lammps import LammpsHook

@pytest.fixture
def simple_atom_lammps_frame():
    h = LammpsHook()
    n_atoms=3
    data_array = h.manipulate_dummy_array("x", n_atoms)
    frame_data = h.lammps_to_frame_data(data_array, positions=True, topology=False)

    return frame_data


def test_topology_lammps_atoms(simple_openmm_topology):
    frame_data = simple_atom_frame(topology=simple_openmm_topology)

    assert len(data.arrays['atom.positions'].index_values.values) == 3


