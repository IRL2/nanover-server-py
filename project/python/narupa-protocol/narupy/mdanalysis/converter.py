import numpy as np
from MDAnalysis import Universe
from MDAnalysis.topology.guessers import guess_atom_element

from narupa.protocol.topology.topology_pb2 import TopologyData
from narupa.protocol.trajectory.frame_pb2 import FrameData

element_index = {
    'H': 1,
    'C': 6,
    'N': 7,
    'O': 8,
    'S': 16
}


def mdanalysis_to_topology_data(u: Universe) -> TopologyData:
    topology_data = TopologyData()

    for residue in u.residues:
        topology_data.arrays['residue.id'].string_values.values.append(residue.resname)
        topology_data.arrays['residue.chain'].index_values.values.append(residue.segment.ix)

    for atom in u.atoms:
        topology_data.arrays['atom.id'].string_values.values.append(atom.name)
        element = element_index[guess_atom_element(atom.name)]
        topology_data.arrays['atom.element'].index_values.values.append(element)
        topology_data.arrays['atom.residue'].index_values.values.append(atom.residue.ix)

    for bond in u.bonds:
        topology_data.arrays['bond'].index_values.values.append(bond.atoms[0].ix)
        topology_data.arrays['bond'].index_values.values.append(bond.atoms[1].ix)

    return topology_data


def mdanalysis_to_frame_data(u: Universe) -> FrameData:
    frame_data = FrameData()

    positions = np.multiply(0.1, np.ndarray.flatten(u.atoms.positions))
    frame_data.arrays["atom.position"].float_values.values.extend(positions)

    return frame_data
