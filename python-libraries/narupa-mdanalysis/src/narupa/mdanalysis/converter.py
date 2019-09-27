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
    """
    Converts from an MDAnalysis universe to Narupa FrameData object.
    :param u: MDAnalysis universe.
    :param topology: Whether to include topology.
    :param positions: Whether to include positions.
    :return: Frame data constructed from MDAnalysis universe.
    """
    frame_data = FrameData()

    if topology:
        frame_data.residue_names = u.residues.resnames
        frame_data.residue_chains = u.residues.segindices
        frame_data.particle_names = u.atoms.names
        frame_data.residue_count = len(u.residues)
        frame_data.chain_count = len(u.segments)
        elements = [element_index[guess_atom_element(name)] for name in u.atoms.names]
        frame_data.particle_elements = elements
        frame_data.particle_residues = u.atoms.resids
        frame_data.bonds = u.atoms.bonds.indices

    if positions:
        frame_data.particle_positions = u.atoms.positions * 0.1

    frame_data.particle_count = len(u.atoms)

    return frame_data
