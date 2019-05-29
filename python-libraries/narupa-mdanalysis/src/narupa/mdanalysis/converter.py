import numpy as np
from MDAnalysis import Universe
from MDAnalysis.topology.guessers import guess_atom_element

from narupa.trajectory import FrameData

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
        frame_data.arrays['residue.id'] = u.residues.resnames
        frame_data.arrays['residue.chain'] = u.residues.segindices
        frame_data.arrays['atom.id'] = u.atoms.names
        elements = [element_index[guess_atom_element(name)] for name in u.atoms.names]
        frame_data.elements = elements
        frame_data.arrays['atom.residue'] = u.atoms.resids
        frame_data.bonds = u.atoms.bonds.indices

    if positions:
        frame_data.particle_positions = u.atoms.positions * 0.1

    return frame_data
