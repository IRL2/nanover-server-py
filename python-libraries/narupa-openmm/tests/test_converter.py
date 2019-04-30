# Pylint does not recognize pytest fixtures, which causes some false warnings.
# pylint: disable=unused-argument,redefined-outer-name
import pytest
from simtk.openmm.app.element import Element
from simtk.openmm.app.topology import Topology

from narupa.openmm import openmm_to_frame_data


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


def test_topology_bonds(simple_openmm_topology):
    data = openmm_to_frame_data(topology=simple_openmm_topology)

    assert len(data.raw.arrays['bond'].index_values.values) == 4


def test_topology_atoms(simple_openmm_topology):
    data = openmm_to_frame_data(topology=simple_openmm_topology)

    assert len(data.elements) == 3
    assert data.elements[0] == 1


def test_topology_residues(simple_openmm_topology):
    data = openmm_to_frame_data(topology=simple_openmm_topology)

    assert len(data.raw.arrays['residue.id'].string_values.values) == 1
