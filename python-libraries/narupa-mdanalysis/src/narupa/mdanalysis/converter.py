import collections

import numpy as np
from MDAnalysis import Universe
from MDAnalysis.topology.guessers import guess_atom_element

from narupa.trajectory import FrameData
from narupa.trajectory.frame_data import PARTICLE_COUNT, RESIDUE_COUNT, CHAIN_COUNT, PARTICLE_ELEMENTS, PARTICLE_NAMES, \
    PARTICLE_RESIDUES, RESIDUE_NAMES, RESIDUE_CHAINS, CHAIN_NAMES, MissingDataError

FrameDataField = collections.namedtuple('FrameDataField', 'key required')

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

FRAME_DATA_TO_MDANALYSIS = {'types': FrameDataField(key=PARTICLE_ELEMENTS, required=True),
                            'names': FrameDataField(PARTICLE_NAMES, False),
                            'resnames': FrameDataField(RESIDUE_NAMES, False),
                            'segnames': FrameDataField(CHAIN_NAMES, False),
                            }

MDA_UNIVERSE_PARAMS_TO_FRAME_DATA = {'n_atoms': PARTICLE_COUNT,
                                     'n_residues': RESIDUE_COUNT,
                                     'n_segments': CHAIN_COUNT,
                                     'atom_resindex': PARTICLE_RESIDUES,
                                     'residue_segindex': RESIDUE_CHAINS}


def mdanalysis_to_frame_data(u: Universe, topology=True, positions=True) -> FrameData:
    """
    Converts from an MDAnalysis universe to Narupa FrameData object.
    :param u: MDAnalysis universe.
    :param topology: Whether to include topology.
    :param positions: Whether to include positions.
    :return: Frame data constructed from MDAnalysis universe.

    :raises: MissingDataError if no positions exist in the MDAnalysis universe, and positions are specified.
    """
    frame_data = FrameData()

    if topology:
        add_mda_attributes(u, frame_data)
        add_mda_counts(u, frame_data)
        add_mda_bonds(u, frame_data)

    if positions:
        add_mda_positions(u, frame_data)

    return frame_data


def frame_data_to_mdanalysis(frame: FrameData):
    """
    Converts from a Narupa FrameData object to an MDAnalysis universe.
    :param u: MDAnalysis universe.
    :param topology: Whether to include topology.
    :param positions: Whether to include positions.
    :return: Frame data constructed from MDAnalysis universe.

    :raises: MissingDataError if no positions exist in the MDAnalysis universe, and positions are specified.
    """
    params = {param_name: try_get_field(frame, field) for param_name, field in
              MDA_UNIVERSE_PARAMS_TO_FRAME_DATA.items()}
    params['trajectory'] = True
    universe = Universe.empty(**params)

    universe.atoms.positions = frame.particle_positions * 10

    add_attributes_to_mda(universe, frame)

    bonds = [(bond[0], bond[1]) for bond in frame.bond_pairs]
    universe.add_TopologyAttr('bonds', bonds)
    return universe


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
    """
    Adds the bonds in a MDAnalysis universe to the frame data, if they exist.
    :param u: MDAnalysis universe.
    :param frame_data: Narupa frame data.
   """
    try:
        frame_data.bonds = u.atoms.bonds.indices
    except AttributeError:
        pass


def add_mda_positions(u: Universe, frame_data: FrameData):
    """
    Adds the positions in a MDAnalysis universe to the frame data, if they exist.
    :param u: MDAnalysis universe.
    :param frame_data: Narupa frame data.

    :raises: MissingDataError, if no positions exist in the universe.
   """
    try:
        frame_data.particle_positions = u.atoms.positions * 0.1
    except AttributeError:
        raise MissingDataError("MDAnalysis universe has no positions.")
    frame_data.particle_count = len(u.atoms)


def try_get_field(frame: FrameData, field):
    return frame.arrays.get(field, None)


def add_attribute_to_mda(universe, attribute_name, values):
    universe.add_TopologyAttr(attribute_name, values)


def add_attributes_to_mda(universe, frame):
    """
    :param universe: MDAnalysis universe to add data to.
    :param frame:
    :raises MissingDataError if a required field is missing.
    :return:
    """
    for name, (key, required) in FRAME_DATA_TO_MDANALYSIS:
        try:
            value = frame.arrays[key]
        except MissingDataError:
            if required:
                raise
            else:
                continue
        universe.add_TopologyAttr(name, value)
