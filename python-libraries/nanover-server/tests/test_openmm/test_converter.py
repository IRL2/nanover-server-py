# Pylint does not recognize pytest fixtures, which causes some false warnings.
# pylint: disable=unused-argument,redefined-outer-name
import itertools

import numpy as np
import pytest
from nanover.openmm import openmm_to_frame_data

from openmm.app.element import Element
from openmm.app.topology import Topology

from simulation_utils import (
    BASIC_SIMULATION_POSITIONS,
    BASIC_SIMULATION_BOX_VECTORS,
    basic_simulation,
)


@pytest.fixture
def simple_openmm_topology():
    topology = Topology()
    chain = topology.addChain(id="A")
    residue = topology.addResidue("RES", chain, "0")
    atom1 = topology.addAtom("Atom1", Element.getByAtomicNumber(1), residue)
    atom2 = topology.addAtom("Atom2", Element.getByAtomicNumber(2), residue)
    atom3 = topology.addAtom("Atom3", Element.getByAtomicNumber(3), residue)
    topology.addBond(atom1, atom2)
    topology.addBond(atom2, atom3)
    return topology


def test_topology_bonds(simple_openmm_topology):
    data = openmm_to_frame_data(topology=simple_openmm_topology)
    assert np.all(data.bond_pairs.flatten() == [0, 1, 1, 2])


def test_topology_atom_elements(simple_openmm_topology):
    data = openmm_to_frame_data(topology=simple_openmm_topology)
    assert np.all(data.particle_elements == [1, 2, 3])


def test_topology_atom_names(simple_openmm_topology):
    data = openmm_to_frame_data(topology=simple_openmm_topology)
    assert data.particle_names == ["Atom1", "Atom2", "Atom3"]


def test_topology_particle_count(simple_openmm_topology):
    data = openmm_to_frame_data(topology=simple_openmm_topology)
    assert data.particle_count == simple_openmm_topology.getNumAtoms()


def test_topology_residues(simple_openmm_topology):
    data = openmm_to_frame_data(topology=simple_openmm_topology)
    assert data.residue_names == ["RES"]


def test_topology_residue_count(simple_openmm_topology):
    data = openmm_to_frame_data(topology=simple_openmm_topology)
    assert data.residue_count == simple_openmm_topology.getNumResidues()


def test_topology_chain_count(simple_openmm_topology):
    data = openmm_to_frame_data(topology=simple_openmm_topology)
    assert data.chain_count == simple_openmm_topology.getNumChains()


def test_topology_chain_names(simple_openmm_topology):
    data = openmm_to_frame_data(topology=simple_openmm_topology)
    assert data.chain_names == ["A"]


def test_topology_particle_residues(simple_openmm_topology):
    data = openmm_to_frame_data(topology=simple_openmm_topology)
    assert np.all(data.particle_residues == [0, 0, 0])


def test_topology_residue_ids(simple_openmm_topology):
    data = openmm_to_frame_data(topology=simple_openmm_topology)
    assert data.residue_ids == ["0"]


def test_topology_residue_chains(simple_openmm_topology):
    data = openmm_to_frame_data(topology=simple_openmm_topology)
    assert np.all(data.residue_chains == [0])


def test_box_vectors(basic_simulation):
    state = basic_simulation.context.getState(getPositions=True)
    data = openmm_to_frame_data(state=state, include_energies=False)
    np.allclose(data.box_vectors, np.asarray(BASIC_SIMULATION_BOX_VECTORS))


def test_positions(basic_simulation):
    state = basic_simulation.context.getState(getPositions=True)
    data = openmm_to_frame_data(state=state, include_energies=False)
    np.allclose(data.particle_positions, np.asarray(BASIC_SIMULATION_POSITIONS))
