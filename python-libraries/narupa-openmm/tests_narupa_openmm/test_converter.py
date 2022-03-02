# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.
# Pylint does not recognize pytest fixtures, which causes some false warnings.
# pylint: disable=unused-argument,redefined-outer-name
import itertools

import numpy as np
import pytest
from numpy.testing import assert_almost_equal
from narupa.openmm import openmm_to_frame_data
from narupa.trajectory import frame_data
try:
    from openmm.app.element import Element
    from openmm.app.topology import Topology
except (ImportError, ModuleNotFoundError):
    from simtk.openmm.app.element import Element
    from simtk.openmm.app.topology import Topology
from .simulation_utils import (
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


# In the following tests, we refer to the raw GRPC FrameData rather than the
# wrapped one to make sure we catch errors due to changes in the wrapper.
# Ultimately, we want to know if we can communicate with the client.
def test_topology_bonds(simple_openmm_topology):
    data = openmm_to_frame_data(topology=simple_openmm_topology)
    assert data.raw.arrays[frame_data.BOND_PAIRS].index_values.values == [0, 1, 1, 2]


def test_topology_atom_elements(simple_openmm_topology):
    data = openmm_to_frame_data(topology=simple_openmm_topology)
    assert data.raw.arrays[frame_data.PARTICLE_ELEMENTS].index_values.values == [1, 2, 3]


def test_topology_atom_names(simple_openmm_topology):
    data = openmm_to_frame_data(topology=simple_openmm_topology)
    assert data.raw.arrays[frame_data.PARTICLE_NAMES].string_values.values == ["Atom1", "Atom2", "Atom3"]


def test_topology_atom_types(simple_openmm_topology):
    """
    Tests that the OpenMM converter does not produce atom types.
    :param simple_openmm_topology: Simple OpenMM topology to test.
    """

    data = openmm_to_frame_data(topology=simple_openmm_topology)
    with pytest.raises(KeyError):
        _ = data.arrays[frame_data.PARTICLE_TYPES].string_values.values


def test_topology_particle_count(simple_openmm_topology):
    data = openmm_to_frame_data(topology=simple_openmm_topology)
    assert data.raw.values[frame_data.PARTICLE_COUNT] == simple_openmm_topology.getNumAtoms()


def test_topology_residues(simple_openmm_topology):
    data = openmm_to_frame_data(topology=simple_openmm_topology)
    assert data.raw.arrays[frame_data.RESIDUE_NAMES].string_values.values == ["RES"]


def test_topology_residue_count(simple_openmm_topology):
    data = openmm_to_frame_data(topology=simple_openmm_topology)
    frame_count = data.raw.values[frame_data.RESIDUE_COUNT]
    frame_count = int(frame_count.number_value)
    assert frame_count == simple_openmm_topology.getNumResidues()


def test_topology_chain_count(simple_openmm_topology):
    data = openmm_to_frame_data(topology=simple_openmm_topology)
    frame_count = data.raw.values[frame_data.CHAIN_COUNT]
    frame_count = int(frame_count.number_value)
    assert frame_count == simple_openmm_topology.getNumChains()


def test_topology_chain_names(simple_openmm_topology):
    data = openmm_to_frame_data(topology=simple_openmm_topology)
    assert data.raw.arrays[frame_data.CHAIN_NAMES].string_values.values == ["A"]


def test_topology_particle_residues(simple_openmm_topology):
    data = openmm_to_frame_data(topology=simple_openmm_topology)
    frame_residues = data.raw.arrays[frame_data.PARTICLE_RESIDUES].index_values.values
    openmm_residues = [0, 0, 0]
    assert frame_residues == openmm_residues


def test_topology_residue_ids(simple_openmm_topology):
    data = openmm_to_frame_data(topology=simple_openmm_topology)
    assert data.raw.arrays[frame_data.RESIDUE_IDS].string_values.values == ["0"]


def test_topology_residue_chains(simple_openmm_topology):
    data = openmm_to_frame_data(topology=simple_openmm_topology)
    frame_chains = data.raw.arrays[frame_data.RESIDUE_CHAINS].index_values.values
    openmm_chains = [0]
    assert frame_chains == openmm_chains


def test_box_vectors(basic_simulation):
    expected = list(itertools.chain(BASIC_SIMULATION_BOX_VECTORS))
    state = basic_simulation.context.getState(getPositions=True)
    data = openmm_to_frame_data(state=state)
    assert_almost_equal(
        data.raw.arrays[frame_data.BOX_VECTORS].float_values.values,
        np.asarray(expected).flatten()
    )


def test_positions(basic_simulation):
    expected = BASIC_SIMULATION_POSITIONS
    state = basic_simulation.context.getState(getPositions=True)
    data = openmm_to_frame_data(state=state)
    assert_almost_equal(
        data.raw.arrays[frame_data.PARTICLE_POSITIONS].float_values.values,
        np.asarray(expected).flatten(),
        decimal=6,
    )


def test_topology_particle_count(simple_openmm_topology):
    data = openmm_to_frame_data(topology=simple_openmm_topology)
    assert data.raw.values[frame_data.PARTICLE_COUNT].number_value == 3
