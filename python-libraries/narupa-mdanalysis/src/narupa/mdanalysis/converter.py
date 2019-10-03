# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.
import collections
from collections import Container

import numpy as np
from MDAnalysis import Universe
from MDAnalysis.topology.guessers import guess_atom_element

from narupa.trajectory import FrameData
from narupa.trajectory.frame_data import (PARTICLE_COUNT, RESIDUE_COUNT, CHAIN_COUNT, PARTICLE_ELEMENTS, PARTICLE_NAMES,
                                          PARTICLE_RESIDUES, RESIDUE_NAMES, RESIDUE_CHAINS, CHAIN_NAMES,
                                          MissingDataError, RESIDUE_IDS)

FrameDataField = collections.namedtuple('FrameDataField', 'key required')
FrameDataFieldConversion = collections.namedtuple('FrameDataFieldConversion', 'key converter')

ELEMENT_INDEX = {
    'H': 1,
    'C': 6,
    'N': 7,
    'O': 8,
    'S': 16,
    'P': 15,
}

INDEX_ELEMENT = {value: key for key, value in ELEMENT_INDEX.items()}

MDANALYSIS_COUNTS_TO_FRAME_DATA = {'atoms': PARTICLE_COUNT,
                                   'residues': RESIDUE_COUNT,
                                   'segments': CHAIN_COUNT,
                                   }
MDANALYSIS_ATOMS_TO_FRAME_DATA = {'types': PARTICLE_ELEMENTS,
                                  'names': PARTICLE_NAMES,
                                  'resindices': PARTICLE_RESIDUES,
                                  }
MDANALYSIS_RESIDUES_TO_FRAME_DATA = {'resnames': RESIDUE_NAMES,
                                     'segindices': RESIDUE_CHAINS,
                                     'resids': RESIDUE_IDS,
                                     }
MDANALYSIS_CHAINS_TO_FRAME_DATA = {'segids': CHAIN_NAMES}

MDANALYSIS_GROUP_TO_ATTRIBUTES = {'atoms': MDANALYSIS_ATOMS_TO_FRAME_DATA,
                                  'residues': MDANALYSIS_RESIDUES_TO_FRAME_DATA,
                                  'segments': MDANALYSIS_CHAINS_TO_FRAME_DATA}

ALL_MDA_ATTRIBUTES = [(group, key, value) for group in MDANALYSIS_GROUP_TO_ATTRIBUTES for key, value in
                      MDANALYSIS_GROUP_TO_ATTRIBUTES[group].items()]


def nullable_int(value):
    if value is None:
        return value
    return int(value)

def _identity(value):
    return value


def _to_chemical_symbol(elements):
    try:
        iterator = iter(elements)
    except TypeError:
        try:
            INDEX_ELEMENT[elements]
        except KeyError:
            raise KeyError(f'Unknown atomic number: {elements}')
    else:
        return [INDEX_ELEMENT[element] for element in elements]


# dictionary of mdanalysis fields to field in frame data, along with any conversion function that needs
# to be applied.
FRAME_DATA_TO_MDANALYSIS = {'types': FrameDataFieldConversion(key=PARTICLE_ELEMENTS, converter=_to_chemical_symbol),
                            'names': FrameDataFieldConversion(PARTICLE_NAMES, _identity),
                            'resnames': FrameDataFieldConversion(RESIDUE_NAMES, _identity),
                            'resids': FrameDataFieldConversion(RESIDUE_IDS, _identity),
                            'segids': FrameDataFieldConversion(CHAIN_NAMES, _identity),
                            }

# dictionary of mdanalysis constructor fields to field in frame data, along with conversion methods.
MDA_UNIVERSE_PARAMS_TO_FRAME_DATA = {'n_atoms': FrameDataFieldConversion(PARTICLE_COUNT, nullable_int),
                                     'n_residues': FrameDataFieldConversion(RESIDUE_COUNT, nullable_int),
                                     'n_segments': FrameDataFieldConversion(CHAIN_COUNT, nullable_int),
                                     'atom_resindex': FrameDataFieldConversion(PARTICLE_RESIDUES, _identity),
                                     'residue_segindex': FrameDataFieldConversion(RESIDUE_CHAINS, _identity)}


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
        _add_mda_attributes(u, frame_data)
        _add_mda_counts_to_frame_data(u, frame_data)
        _add_mda_bonds_to_frame_data(u, frame_data)

    if positions:
        _add_mda_positions_to_frame_data(u, frame_data)

    return frame_data


def frame_data_to_mdanalysis(frame: FrameData) -> Universe:
    """
    Converts from a Narupa FrameData object to an MDAnalysis universe.

    :param frame: Narupa FrameData object.
    :return: MDAnalysis universe constructed from the given FrameData.
    """

    params = _get_universe_constructor_params(frame)
    universe = Universe.empty(**params)

    add_positions_to_mda(universe, frame)

    # additional topology information.
    _add_frame_attributes_to_mda(universe, frame)
    add_bonds_to_mda(universe, frame)

    return universe


def add_positions_to_mda(u: Universe, frame: FrameData):
    u.atoms.positions = np.array(frame.particle_positions) * 10


def add_bonds_to_mda(u: Universe, frame: FrameData):
    """
    Add bonds from a framedata object to an MDAnalysis universe.
    :param u: MDAnalysis universe.
    :param frame: Narupa FrameData.
    """
    try:
        bonds = [(bond[0], bond[1]) for bond in frame.bonds]
    except MissingDataError:
        return
    u.add_TopologyAttr('bonds', bonds)


def _get_universe_constructor_params(frame: FrameData):
    """
    Gets the MDAnalysis universe constructor params from a Narupa frame data.
    :param frame: Narupa FrameData object.
    :return: Dictionary of params to construct an MDAnalysis universe object.

    The MDAnalysis universe empty constructor takes several optional parameters used to define
    options such as number of atoms, number of residues, number of segments, and their identifiers.
    This method extracts this data from a Narupa FrameData object.
    """
    params = {
        param_name: converter(_try_get_field(frame, field))
        for param_name, (field, converter)
        in MDA_UNIVERSE_PARAMS_TO_FRAME_DATA.items()
    }

    params['trajectory'] = True
    return params


def get_mda_attribute(u: Universe, group, group_attribute):
    """
    Gets an attribute associated with a particular group.
    :param u: MDAnalysis universe.
    :param group: The group in the MDAnalysis universe in which the attribute exists.
    :param group_attribute: The attribute.
    :return: The attribute, if it exists.
    :raises: AttributeError: If either the universe does not contain the given group, or the attribute
    does not exist in the given group, an AttributeError will be raised.
    """
    return getattr(getattr(u, group), group_attribute)


def _add_mda_attributes(u: Universe, frame_data: FrameData):
    """
    Adds all available MDAnalysis attributes from the given universe to the given frame data
    :param u: MDAnalysis universe.
    :param frame_data: Narupa frame data.

    Adds particle, residue and chain information, if available.
    """
    for group, attribute, frame_key in ALL_MDA_ATTRIBUTES:
        try:
            field = get_mda_attribute(u, group, attribute)
        except AttributeError:
            continue
        if frame_key == PARTICLE_ELEMENTS:
            field = [ELEMENT_INDEX[guess_atom_element(name)] for name in field]
        frame_data.arrays[frame_key] = field


def _add_mda_counts_to_frame_data(u: Universe, frame_data: FrameData):
    """
    Adds the counts of all available MDAnalysis groups from the given universe to the given frame data.
    :param u: MDAnalysis universe.
    :param frame_data: Narupa frame data.

    Adds particle counts, residue counts and chain counts, if available.
    """
    for attribute, frame_key in MDANALYSIS_COUNTS_TO_FRAME_DATA.items():
        try:
            field = getattr(u, attribute)
        except AttributeError:
            continue
        frame_data.values[frame_key] = len(field)


def _add_mda_bonds_to_frame_data(u: Universe, frame_data: FrameData):
    """
    Adds the bonds in a MDAnalysis universe to the frame data, if they exist.
    :param u: MDAnalysis universe.
    :param frame_data: Narupa frame data.
   """
    try:
        frame_data.bonds = u.atoms.bonds.indices
    except AttributeError:
        pass


def _add_mda_positions_to_frame_data(u: Universe, frame_data: FrameData):
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


def _try_get_field(frame: FrameData, field):
    array_keys = frame.array_keys
    value_keys = frame.value_keys
    if field in array_keys:
        return frame.arrays.get(field)
    elif field in value_keys:
        return frame.values.get(field)


def _add_frame_attributes_to_mda(universe, frame):
    for name, (key, converter) in FRAME_DATA_TO_MDANALYSIS.items():
        try:
            value = frame.arrays[key]
        except (KeyError, MissingDataError):
            # TODO should some fields be required?
            continue
        universe.add_TopologyAttr(name, converter(value))
