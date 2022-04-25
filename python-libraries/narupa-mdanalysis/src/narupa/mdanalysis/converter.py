# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.
"""
Module for performing conversions between MDAnalysis universes and Narupa FrameData objects.
"""
import collections

import numpy as np
from MDAnalysis import Universe
from MDAnalysis.topology.guessers import guess_atom_element

from narupa.trajectory import FrameData
from narupa.trajectory.frame_data import (
    PARTICLE_COUNT, RESIDUE_COUNT, CHAIN_COUNT,
    PARTICLE_ELEMENTS, PARTICLE_NAMES, PARTICLE_RESIDUES,
    RESIDUE_NAMES, RESIDUE_CHAINS, RESIDUE_IDS,
    CHAIN_NAMES,
    MissingDataError,
)

# tuple for storing a frame data key and whether it is required in conversion.
FrameDataField = collections.namedtuple('FrameDataField', 'key required')
# tuple for storing a frame data key and a conversion method to apply when
# producing the corresponding attribute in MDAnalysis.
FrameDataFieldConversion = collections.namedtuple('FrameDataFieldConversion', 'key converter')

ELEMENT_INDEX = {
    'H': 1,
    'He': 2,
    'Li': 3,
    'Be': 4,
    'B': 5,
    'C': 6,
    'N': 7,
    'O': 8,
    'F': 9,
    'Ne': 10,
    'Na': 11,
    'Mg': 12,
    'Al': 13,
    'Si': 14,
    'P': 15,
    'S': 16,
    'Cl': 17,
    'Ar': 18,
    'K': 19,
    'Ca': 20,
    'Sc': 21,
    'Ti': 22,
    'V': 23,
    'Cr': 24,
    'Mn': 25,
    'Fe': 26,
    'Co': 27,
    'Ni': 28,
    'Cu': 29,
    'Zn': 30,
    'Ga': 31,
    'Ge': 32,
    'As': 33,
    'Se': 34,
    'Br': 35,
    'Kr': 36,
    'Rb': 37,
    'Sr': 38,
    'Y': 39,
    'Zr': 40,
    'Nb': 41,
    'Mo': 42,
    'Tc': 43,
    'Ru': 44,
    'Rh': 45,
    'Pd': 46,
    'Ag': 47,
    'Cd': 48,
    'In': 49,
    'Sn': 50,
    'Sb': 51,
    'Te': 52,
    'I': 53,
    'Xe': 54,
    'Cs': 55,
    'Ba': 56,
    'La': 57,
    'Ce': 58,
    'Pr': 59,
    'Nd': 60,
    'Pm': 61,
    'Sm': 62,
    'Eu': 63,
    'Gd': 64,
    'Tb': 65,
    'Dy': 66,
    'Ho': 67,
    'Er': 68,
    'Tu': 69,
    'Yb': 70,
    'Lu': 71,
    'Hf': 72,
    'Ta': 73,
    'W': 74,
    'Re': 75,
    'Os': 76,
    'Ir': 77,
    'Pt': 78,
    'Au': 79,
    'Hg': 80,
    'Tl': 81,
    'Pb': 82,
    'Bi': 83,
    'Po': 84,
    'At': 85,
    'Rn': 86,
    'Fr': 87,
    'Ra': 88,
    'Ac': 89,
    'Th': 90,
    'Pa': 91,
    'U': 92,
    'Np': 93,
    'Pu': 94,
    'Am': 95,
    'Cm': 96,
    'Bk': 97,
    'Cf': 98,
    'Es': 99,
    'Fm': 100,
    'Md': 101,
    'No': 102,
    'Lr': 103,
    'Rf': 104,
    'Db': 105,
    'Sg': 106,
    'Bh': 107,
    'Hs': 108,
    'Mt': 109,
    'Ds': 110,
    'Rg': 111,
    'Cn': 112,
    'Nh': 113,
    'Fv': 114,
    'Ms': 115,
    'Lv': 116,
    'Ts': 117,
    'Og': 118
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

ALL_MDA_ATTRIBUTES = [
    (group, key, value)
    for group in MDANALYSIS_GROUP_TO_ATTRIBUTES
    for key, value in MDANALYSIS_GROUP_TO_ATTRIBUTES[group].items()
]


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

    :param u: MDAnalysis :class:`Universe`.
    :param topology: Whether to include topology.
    :param positions: Whether to include positions.
    :return: :class:`FrameData` constructed from MDAnalysis universe.

    :raises MissingDataError: if no positions exist in the MDAnalysis universe,
        and positions are specified.

    Topological information consists any available information such as bonds,
    residue names, residue ids, atom names, chain names, residue index and
    chain indexes
    """
    frame_data = FrameData()

    if topology:
        add_mda_topology_to_frame_data(u, frame_data)

    if positions:
        add_mda_positions_to_frame_data(u, frame_data)

    return frame_data


def frame_data_to_mdanalysis(frame: FrameData) -> Universe:
    """
    Converts from a Narupa :class:`FrameData` object to an MDAnalysis universe.

    :param frame: Narupa :class:`FrameData` object.
    :return: MDAnalysis :class:`Universe` constructed from the given FrameData.
    """

    params = _get_universe_constructor_params(frame)
    universe = Universe.empty(**params)

    add_frame_positions_to_mda(universe, frame)

    # additional topology information.
    _add_frame_attributes_to_mda(universe, frame)
    _add_bonds_to_mda(universe, frame)

    return universe


def add_mda_topology_to_frame_data(u, frame_data):
    """
    Adds available topology information from an MDAnalysis Universe to a FrameData.
    
    :param u: MDAnalysis :class:`Universe`.
    :param frame_data: :class:`FrameData` to add to.
    """
    _add_mda_attributes_to_frame_data(u, frame_data)
    _add_mda_counts_to_frame_data(u, frame_data)
    _add_mda_bonds_to_frame_data(u, frame_data)


def add_mda_positions_to_frame_data(u: Universe, frame_data: FrameData):
    """
    Adds the positions in a MDAnalysis universe to the frame data, if they exist.
    
    :param u: MDAnalysis :class:`Universe`.
    :param frame_data: Narupa :class:`FrameData` to add to.

    :raises MissingDataError: if no positions exist in the universe.
   """
    try:
        frame_data.particle_positions = u.atoms.positions * 0.1
    except AttributeError:
        raise MissingDataError("MDAnalysis universe has no positions.")


def add_frame_topology_to_mda(u: Universe, frame: FrameData):
    _add_bonds_to_mda(u, frame)
    _add_frame_attributes_to_mda(u, frame)


def add_frame_positions_to_mda(u: Universe, frame: FrameData):
    """
    Updates the positions in an MDAnalysis :class:`Universe` with those from the given frame.
    
    :param u: MDAnalysis :class:`Universe` to set positions of.
    :param frame: Narupa :class:`FrameData` from which to extract positions.
    """
    u.atoms.positions = np.array(frame.particle_positions) * 10


def _add_bonds_to_mda(u: Universe, frame: FrameData):
    """
    Add bonds from a framedata object to an MDAnalysis universe.

    :param u: MDAnalysis :class:`Universe`.
    :param frame: Narupa :class:`FrameData`.
    """
    try:
        bonds = [(bond[0], bond[1]) for bond in frame.bond_pairs]
    except MissingDataError:
        return
    u.add_TopologyAttr('bonds', bonds)


def _get_universe_constructor_params(frame: FrameData):
    """
    Gets the MDAnalysis universe constructor params from a Narupa frame data.

    :param frame: Narupa FrameData object.
    :return: Dictionary of params to construct an MDAnalysis universe object.

    The MDAnalysis universe empty constructor takes several optional parameters
    used to define options such as number of atoms, number of residues, number
    of segments, and their identifiers. This method extracts this data from a
    Narupa :class:`FrameData` object.
    """
    params = {
        param_name: converter(_try_get_field(frame, field))
        for param_name, (field, converter)
        in MDA_UNIVERSE_PARAMS_TO_FRAME_DATA.items()
    }

    # strip unused arguments
    params = {key: value for key, value in params.items() if value is not None}
    params['trajectory'] = True

    if 'atom_resindex' not in params and 'n_atoms' in params:
        params['atom_resindex'] = [0] * params['n_atoms']
    if 'residue_segindex' not in params and 'atom_resindex' in params:
        n_residues = params.get('n_residues', max(params['atom_resindex']) + 1)
        params['residue_segindex'] = [0] * n_residues

    return params


def _get_mda_attribute(u: Universe, group, group_attribute):
    """
    Gets an attribute associated with a particular group.

    :param u: MDAnalysis universe.
    :param group: The group in the MDAnalysis universe in which the attribute exists.
    :param group_attribute: The attribute.
    :return: The attribute, if it exists.
    :raises AttributeError: If either the universe does not contain the given
        group, or the attribute does not exist in the given group, an
        :exc:`AttributeError` will be raised.
    """
    return getattr(getattr(u, group), group_attribute)


def _add_mda_attributes_to_frame_data(u: Universe, frame_data: FrameData):
    """
    Adds all available MDAnalysis attributes from the given universe to the given frame data

    :param u: MDAnalysis universe.
    :param frame_data: Narupa frame data.

    Adds particle, residue and chain information, if available.
    """
    for group, attribute, frame_key in ALL_MDA_ATTRIBUTES:
        try:
            field = _get_mda_attribute(u, group, attribute)
        except AttributeError:
            continue
        if frame_key == PARTICLE_ELEMENTS:
            # When MDAnalysis guesses an element symbol, it returns it fully
            # in upper case. We need to fix the case before we can query our
            # table.
            field = [
                ELEMENT_INDEX[guess_atom_element(name).capitalize()]
                for name in field
            ]
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
        frame_data.bond_pairs = u.atoms.bonds.indices
    except AttributeError:
        pass


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
