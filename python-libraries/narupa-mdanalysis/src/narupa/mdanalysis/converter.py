import numpy as np
from MDAnalysis import Universe
from MDAnalysis.topology.guessers import guess_atom_element

from narupa.trajectory import FrameData
from narupa.trajectory.frame_data import PARTICLE_COUNT, RESIDUE_COUNT, CHAIN_COUNT, PARTICLE_ELEMENTS, PARTICLE_NAMES, \
    PARTICLE_RESIDUES, RESIDUE_NAMES, RESIDUE_CHAINS, CHAIN_NAMES, MissingDataError

ELEMENT_INDEX = {
    'H': 1,
    'C': 6,
    'N': 7,
    'O': 8,
    'S': 16,
    'P': 15,
}

INDEX_ELEMENT = {value: key for key, value in ELEMENT_INDEX.items()}

FRAME_DATA_TO_MDANALYSIS_COUNTS = {'atoms': PARTICLE_COUNT,
                                   'residues': RESIDUE_COUNT,
                                   'segments': CHAIN_COUNT,
                                   }
FRAME_DATA_TO_MDANALYSIS_ATOMS = {'types': PARTICLE_ELEMENTS,
                                  'names': PARTICLE_NAMES,
                                  'resids': PARTICLE_RESIDUES,
                                  }
FRAME_DATA_TO_MDANALYSIS_RESIDUES = {'resnames': RESIDUE_NAMES,
                                     'segindices': RESIDUE_CHAINS,
                                     }
FRAME_DATA_TO_MDANALYSIS_CHAINS = {'segids': CHAIN_NAMES}
GROUP_TO_ATTRIBUTE = {'atoms': FRAME_DATA_TO_MDANALYSIS_ATOMS, 'residues': FRAME_DATA_TO_MDANALYSIS_RESIDUES,
                      'segments': FRAME_DATA_TO_MDANALYSIS_CHAINS}
ALL_MDA_ATTRIBUTES = [(group, key, value) for group in GROUP_TO_ATTRIBUTE for key, value in
                      GROUP_TO_ATTRIBUTE[group].items()]


def get_attribute(u: Universe, universe_attribute, group_attribute):
    return getattr(getattr(u, universe_attribute), group_attribute)


def add_mda_attributes(u: Universe, frame_data: FrameData):
    """
    Adds all available MDAnalysis attributes from the given universe to the given frame data
    :param u: MDAnalysis universe.
    :param frame_data: Narupa frame data.

    Adds particle, residue and chain information, if available.
    """
    for group, attribute, frame_key in ALL_MDA_ATTRIBUTES:
        try:
            field = get_attribute(u, group, attribute)
        except AttributeError:
            continue
        if frame_key == PARTICLE_ELEMENTS:
            field = [ELEMENT_INDEX[guess_atom_element(name)] for name in field]
        frame_data.arrays[frame_key] = field


def add_mda_counts(u: Universe, frame_data: FrameData):
    """
    Adds the counts of all available MDAnalysis groups from the given universe to the given frame data.
    :param u: MDAnalysis universe.
    :param frame_data: Narupa frame data.

    Adds particle counts, residue counts and chain counts, if available.
    """
    for attribute, frame_key in FRAME_DATA_TO_MDANALYSIS_COUNTS.items():
        try:
            field = getattr(u, attribute)
        except AttributeError:
            continue
        frame_data.values[frame_key] = len(field)


def add_mda_bonds(u: Universe, frame_data: FrameData):
    try:
        frame_data.bonds = u.atoms.bonds.indices
    except AttributeError:
        pass


def add_mda_positions(u: Universe, frame_data: FrameData):
    try:
        frame_data.particle_positions = u.atoms.positions * 0.1
    except AttributeError:
        raise MissingDataError("MDAnalysis universe has no positions.")
    frame_data.particle_count = len(u.atoms)


def mdanalysis_to_frame_data(u: Universe, topology=True, positions=True) -> FrameData:
    """
    Converts from an MDAnalysis universe to Narupa FrameData object.
    :param u: MDAnalysis universe.
    :param topology: Whether to include topology.
    :param positions: Whether to include positions.
    :return: Frame data constructed from MDAnalysis universe.

    :raises: AttributeError
    """
    frame_data = FrameData()

    if topology:
        add_mda_attributes(u, frame_data)
        add_mda_counts(u, frame_data)
        add_mda_bonds(u, frame_data)

    if positions:
        add_mda_positions(u, frame_data)

    return frame_data
