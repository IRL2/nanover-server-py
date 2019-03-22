from narupa.protocol.topology.topology_pb2 import TopologyData

from simtk.openmm.app.topology import Topology


def openmm_topology_to_topology_data(topology : Topology) -> TopologyData:
    data = TopologyData()

    data.arrays['residue.id'].string_values.values.extend([residue.name for residue in topology.residues()])
    data.arrays['residue.chain'].index_values.values.extend([residue.chain.index for residue in topology.residues()])

    for atom in topology.atoms():
        data.arrays['atom.id'].string_values.values.append(atom.name)
        data.arrays['atom.element'].index_values.values.append(atom.element.atomic_number)
        data.arrays['atom.residue'].index_values.values.append(atom.residue.index)

    for bond in topology.bonds():
        data.arrays['bond'].index_values.values.append(bond[0].index)
        data.arrays['bond'].index_values.values.append(bond[1].index)

    return data