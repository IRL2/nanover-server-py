"""
Module for performing conversions between MDAnalysis universes and NanoVer FrameData objects.
"""

import collections
from contextlib import suppress

import numpy as np
from MDAnalysis import Universe
from MDAnalysis.guesser.default_guesser import DefaultGuesser

from nanover.trajectory import FrameData, FrameData2
from nanover.trajectory.frame_data import (
    PARTICLE_COUNT,
    RESIDUE_COUNT,
    CHAIN_COUNT,
    PARTICLE_ELEMENTS,
    PARTICLE_NAMES,
    PARTICLE_RESIDUES,
    RESIDUE_NAMES,
    RESIDUE_CHAINS,
    RESIDUE_IDS,
    CHAIN_NAMES,
    MissingDataError,
)
from nanover.trajectory.convert import convert_dict_frame_to_grpc_frame

# tuple for storing a frame data key and whether it is required in conversion.
FrameDataField = collections.namedtuple("FrameDataField", "key required")
# tuple for storing a frame data key and a conversion method to apply when
# producing the corresponding attribute in MDAnalysis.
FrameDataFieldConversion = collections.namedtuple(
    "FrameDataFieldConversion", "key converter"
)

ELEMENT_NAMES = ",H,He,Li,Be,B,C,N,O,F,Ne,Na,Mg,Al,Si,P,S,Cl,Ar,K,Ca,Sc,Ti,V,Cr,Mn,Fe,Co,Ni,Cu,Zn,Ga,Ge,As,Se,Br,Kr,Rb,Sr,Y,Zr,Nb,Mo,Tc,Ru,Rh,Pd,Ag,Cd,In,Sn,Sb,Te,I,Xe,Cs,Ba,La,Ce,Pr,Nd,Pm,Sm,Eu,Gd,Tb,Dy,Ho,Er,Tu,Yb,Lu,Hf,Ta,W,Re,Os,Ir,Pt,Au,Hg,Tl,Pb,Bi,Po,At,Rn,Fr,Ra,Ac,Th,Pa,U,Np,Pu,Am,Cm,Bk,Cf,Es,Fm,Md,No,Lr,Rf,Db,Sg,Bh,Hs,Mt,Ds,Rg,Cn,Nh,Fv,Ms,Lv,Ts,Og"
ELEMENT_INDEX = {name: index for index, name in enumerate(ELEMENT_NAMES.split(","))}
INDEX_ELEMENT = {index: name for name, index in ELEMENT_INDEX.items()}

MDANALYSIS_COUNTS_TO_FRAME_DATA = {
    "atoms": PARTICLE_COUNT,
    "residues": RESIDUE_COUNT,
    "segments": CHAIN_COUNT,
}
MDANALYSIS_ATOMS_TO_FRAME_DATA = {
    "types": PARTICLE_ELEMENTS,
    "names": PARTICLE_NAMES,
    "resindices": PARTICLE_RESIDUES,
}
MDANALYSIS_RESIDUES_TO_FRAME_DATA = {
    "resnames": RESIDUE_NAMES,
    "segindices": RESIDUE_CHAINS,
    "resids": RESIDUE_IDS,
}
MDANALYSIS_CHAINS_TO_FRAME_DATA = {"segids": CHAIN_NAMES}

MDANALYSIS_GROUP_TO_ATTRIBUTES = {
    "atoms": MDANALYSIS_ATOMS_TO_FRAME_DATA,
    "residues": MDANALYSIS_RESIDUES_TO_FRAME_DATA,
    "segments": MDANALYSIS_CHAINS_TO_FRAME_DATA,
}

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
        iter(elements)
    except TypeError:
        try:
            INDEX_ELEMENT[elements]
        except KeyError:
            raise KeyError(f"Unknown atomic number: {elements}")
    else:
        return [INDEX_ELEMENT[element] for element in elements]


# dictionary of mdanalysis fields to field in frame data, along with any conversion function that needs
# to be applied.
FRAME_DATA_TO_MDANALYSIS = {
    "types": FrameDataFieldConversion(
        key=PARTICLE_ELEMENTS, converter=_to_chemical_symbol
    ),
    "names": FrameDataFieldConversion(PARTICLE_NAMES, _identity),
    "resnames": FrameDataFieldConversion(RESIDUE_NAMES, _identity),
    "resids": FrameDataFieldConversion(RESIDUE_IDS, _identity),
    "segids": FrameDataFieldConversion(CHAIN_NAMES, _identity),
}

# dictionary of mdanalysis constructor fields to field in frame data, along with conversion methods.
MDA_UNIVERSE_PARAMS_TO_FRAME_DATA = {
    "n_atoms": FrameDataFieldConversion(PARTICLE_COUNT, nullable_int),
    "n_residues": FrameDataFieldConversion(RESIDUE_COUNT, nullable_int),
    "n_segments": FrameDataFieldConversion(CHAIN_COUNT, nullable_int),
    "atom_resindex": FrameDataFieldConversion(PARTICLE_RESIDUES, _identity),
    "residue_segindex": FrameDataFieldConversion(RESIDUE_CHAINS, _identity),
}


def mdanalysis_to_frame_data(u: Universe, topology=True, positions=True) -> FrameData:
    """
    Converts from an MDAnalysis universe to NanoVer FrameData object.

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


def frame_data_to_mdanalysis(frame: FrameData | FrameData2) -> Universe:
    """
    Converts from a NanoVer :class:`FrameData` object to an MDAnalysis universe.

    :param frame: NanoVer :class:`FrameData` object.
    :return: MDAnalysis :class:`Universe` constructed from the given FrameData.
    """
    if isinstance(frame, FrameData2):
        frame = convert_dict_frame_to_grpc_frame(frame.frame_dict)

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
    :param frame_data: NanoVer :class:`FrameData` to add to.

    :raises MissingDataError: if no positions exist in the universe.
    """
    try:
        # convert from angstroms (mdanalysis) to nanometes (nanover)
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
    :param frame: NanoVer :class:`FrameData` from which to extract positions.
    """
    # convert from nanometers (nanover) to angstroms (mdanalysis)
    u.atoms.positions = np.array(frame.particle_positions) * 10


def _add_bonds_to_mda(u: Universe, frame: FrameData):
    """
    Add bonds from a framedata object to an MDAnalysis universe.

    :param u: MDAnalysis :class:`Universe`.
    :param frame: NanoVer :class:`FrameData`.
    """
    with suppress(MissingDataError):
        # TODO: why does mypy hate this?
        bonds = [(bond[0], bond[1]) for bond in frame.bond_pairs]  # type: ignore
        u.add_TopologyAttr("bonds", bonds)


def _get_universe_constructor_params(frame: FrameData):
    """
    Gets the MDAnalysis universe constructor params from a NanoVer frame data.

    :param frame: NanoVer FrameData object.
    :return: Dictionary of params to construct an MDAnalysis universe object.

    The MDAnalysis universe empty constructor takes several optional parameters
    used to define options such as number of atoms, number of residues, number
    of segments, and their identifiers. This method extracts this data from a
    NanoVer :class:`FrameData` object.
    """
    params = {
        param_name: converter(_try_get_field(frame, field))
        for param_name, (field, converter) in MDA_UNIVERSE_PARAMS_TO_FRAME_DATA.items()
    }

    # strip unused arguments
    params = {key: value for key, value in params.items() if value is not None}
    params["trajectory"] = True

    if "atom_resindex" not in params and "n_atoms" in params:
        params["atom_resindex"] = [0] * params["n_atoms"]
    if "residue_segindex" not in params and "atom_resindex" in params:
        n_residues = params.get("n_residues", max(params["atom_resindex"]) + 1)
        params["residue_segindex"] = [0] * n_residues

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
    :param frame_data: NanoVer frame data.

    Adds particle, residue and chain information, if available.
    """
    guesser = DefaultGuesser(u)

    for group, attribute, frame_key in ALL_MDA_ATTRIBUTES:
        with suppress(AttributeError):
            field = _get_mda_attribute(u, group, attribute)
            if frame_key == PARTICLE_ELEMENTS:
                # When MDAnalysis guesses an element symbol, it returns it fully
                # in upper case. We need to fix the case before we can query our
                # table.
                field = [
                    ELEMENT_INDEX[guesser.guess_atom_element(name).capitalize()]
                    for name in field
                ]
            frame_data.arrays[frame_key] = field


def _add_mda_counts_to_frame_data(u: Universe, frame_data: FrameData):
    """
    Adds the counts of all available MDAnalysis groups from the given universe to the given frame data.

    :param u: MDAnalysis universe.
    :param frame_data: NanoVer frame data.

    Adds particle counts, residue counts and chain counts, if available.
    """
    for attribute, frame_key in MDANALYSIS_COUNTS_TO_FRAME_DATA.items():
        with suppress(AttributeError):
            field = getattr(u, attribute)
            frame_data.values[frame_key] = len(field)


def _add_mda_bonds_to_frame_data(u: Universe, frame_data: FrameData):
    """
    Adds the bonds in a MDAnalysis universe to the frame data, if they exist.

    :param u: MDAnalysis universe.
    :param frame_data: NanoVer frame data.
    """
    with suppress(AttributeError):
        frame_data.bond_pairs = u.atoms.bonds.indices


def _try_get_field(frame: FrameData, field):
    array_keys = frame.array_keys
    value_keys = frame.value_keys
    if field in array_keys:
        return frame.arrays.get(field)
    elif field in value_keys:
        return frame.values.get(field)


def _add_frame_attributes_to_mda(universe, frame):
    for name, (key, converter) in FRAME_DATA_TO_MDANALYSIS.items():
        with suppress(KeyError, MissingDataError):
            value = frame.arrays[key]
            universe.add_TopologyAttr(name, converter(value))
