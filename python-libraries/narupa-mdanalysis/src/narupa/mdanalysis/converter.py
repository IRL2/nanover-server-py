import numpy as np
from MDAnalysis import Universe
from MDAnalysis.topology.guessers import guess_atom_element

from narupa.protocol.trajectory import FrameData

element_index = {
    'H': 1,
    'C': 6,
    'N': 7,
    'O': 8,
    'S': 16
}


def mdanalysis_to_frame_data(u: Universe, topology=True, positions=True) -> FrameData:
    frame_data = FrameData()

    if topology:
        for residue in u.residues:
            frame_data.arrays['residue.id'].string_values.values.append(residue.resname)
            frame_data.arrays['residue.chain'].index_values.values.append(residue.segment.ix)

        for atom in u.atoms:
            frame_data.arrays['atom.id'].string_values.values.append(atom.name)
            element = element_index[guess_atom_element(atom.name)]
            frame_data.arrays['atom.element'].index_values.values.append(element)
            frame_data.arrays['atom.residue'].index_values.values.append(atom.residue.ix)

        for bond in u.bonds:
            frame_data.arrays['bond'].index_values.values.append(bond.atoms[0].ix)
            frame_data.arrays['bond'].index_values.values.append(bond.atoms[1].ix)

    if positions:
        positions = np.multiply(0.1, np.ndarray.flatten(u.atoms.positions))
        frame_data.arrays["atom.position"].float_values.values.extend(positions)

    return frame_data