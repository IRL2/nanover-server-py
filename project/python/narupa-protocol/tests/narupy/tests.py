"""
Copyright (c) Alex Jamieson-Binnie, University Of Bristol. All rights reserved.
Licensed under the GPL. See License.txt in the project root for license information.
"""

import pytest
from simtk.openmm.app.element import Element
from simtk.openmm.app.topology import Topology

from narupy.openmm import openmm_to_frame_data


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

    assert len(data.arrays['bond'].index_values.values) == 4
