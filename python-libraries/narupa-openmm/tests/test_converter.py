# Pylint does not recognize pytest fixtures, which causes some false warnings.
# pylint: disable=unused-argument,redefined-outer-name
import pytest
from simtk.openmm.app.element import Element
from simtk.openmm.app.topology import Topology

from narupa.openmm import openmm_to_frame_data
from narupa.trajectory import frame_data


@pytest.fixture
def simple_openmm_topology():
    topology = Topology()
    chain = topology.addChain()
    residue = topology.addResidue("RES", chain, 0)
    atom1 = topology.addAtom("Atom1", Element.getByAtomicNumber(1), residue)
    atom2 = topology.addAtom("Atom2", Element.getByAtomicNumber(2), residue)
    atom3 = topology.addAtom("Atom3", Element.getByAtomicNumber(3), residue)
    topology.addBond(atom1, atom2)
    topology.addBond(atom2, atom3)
    return topology


# In the following tests, we refer to the raw GRPC FrameData rather than the
# wrapped one to make sure we catch errors due to changes in the wrapper.
# Ultimately, we want to know if we can communicate with the client.
def test_topology_bonds(simple_openmm_topology):
    data = openmm_to_frame_data(topology=simple_openmm_topology)

    assert data.raw.arrays[frame_data.BOND_PAIRS].index_values.values == [0, 1, 1, 2]


def test_topology_atoms(simple_openmm_topology):
    data = openmm_to_frame_data(topology=simple_openmm_topology)

    assert data.raw.arrays[frame_data.PARTICLE_ELEMENTS].index_values.values == [1, 2, 3]


def test_topology_residues(simple_openmm_topology):
    data = openmm_to_frame_data(topology=simple_openmm_topology)

    assert data.raw.arrays[frame_data.RESIDUE_NAMES].string_values.values == ["RES"]


def test_topology_particle_count(simple_openmm_topology):
    data = openmm_to_frame_data(topology=simple_openmm_topology)

    assert data.raw.values[frame_data.PARTICLE_COUNT].number_value == 3