# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.
import os

import numpy as np
import pytest
from MDAnalysis import Universe
from narupa.mdanalysis.converter import (
    INDEX_ELEMENT,
    MDANALYSIS_COUNTS_TO_FRAME_DATA,
    ALL_MDA_ATTRIBUTES,
    mdanalysis_to_frame_data,
    frame_data_to_mdanalysis,
    add_mda_topology_to_frame_data,
    _get_mda_attribute,
)
from narupa.trajectory.frame_data import (PARTICLE_ELEMENTS, MissingDataError, FrameData)

TEST_SYSTEM_PATH = os.path.join(
    os.path.dirname(os.path.realpath(__file__)),
    '2efv_fragment.pdb',
)


@pytest.fixture
def universe():
    return Universe(TEST_SYSTEM_PATH, guess_bonds=True)


@pytest.fixture()
def empty_universe():
    return Universe.empty(
        n_atoms=1,
        atom_resindex=[0],
        residue_segindex=[0],
        trajectory=True,
    )


@pytest.fixture()
def empty_universe_no_positions():
    return Universe.empty(
        n_atoms=1,
        atom_resindex=[0],
        residue_segindex=[0],
        trajectory=False,
    )


@pytest.fixture()
def single_atom_universe(empty_universe: Universe):
    empty_universe.atoms.positions = [[0, 0, 0]]
    return empty_universe


@pytest.fixture()
def metal_universe():
    """
    This universe only has atom types. These types correspond to element symbols
    that have more than one letter.
    """
    universe = Universe.empty(
        n_atoms=3,
        atom_resindex=[0, 0, 0],
        residue_segindex=[0],
    )
    universe.add_TopologyAttr('types', ['NA', 'CL', 'FE'])
    return universe


@pytest.fixture
def frame_data_and_universe(universe):
    return mdanalysis_to_frame_data(universe), universe


def test_mdanalysis_to_frame_data(universe):
    frame_data = mdanalysis_to_frame_data(universe)
    assert frame_data is not None


@pytest.mark.parametrize(
    "universe_attribute, mda_attribute, frame_field",
    ALL_MDA_ATTRIBUTES,
)
def test_mdanalysis_particle_field(universe_attribute, mda_attribute, frame_field, frame_data_and_universe):
    frame, universe = frame_data_and_universe
    # fetches the atoms, residues or chains object, then the attribute.
    attribute = _get_mda_attribute(universe, universe_attribute, mda_attribute)
    field = frame.arrays[frame_field]
    if frame_field == PARTICLE_ELEMENTS:
        field = [INDEX_ELEMENT[x] for x in field]
    assert all(a == b for a, b in zip(attribute, field))


def test_mdanalysis_positions(frame_data_and_universe):
    frame, universe = frame_data_and_universe
    assert np.allclose(np.array(frame.particle_positions) * 10, universe.atoms.positions)


@pytest.mark.parametrize(
    "mda_attribute, frame_field",
    list(MDANALYSIS_COUNTS_TO_FRAME_DATA.items())
)
def test_mdanalysis_counts(mda_attribute, frame_field, frame_data_and_universe):
    frame, universe = frame_data_and_universe
    assert len(getattr(universe, mda_attribute)) == frame.values[frame_field]


def test_mdanalysis_bonds(frame_data_and_universe):
    frame, universe = frame_data_and_universe
    frame_bonds = np.array(frame.bond_pairs)
    universe_bonds = np.array(universe.bonds.indices)
    assert np.all(frame_bonds == universe_bonds)


def test_empty_universe(empty_universe_no_positions):
    with pytest.raises(MissingDataError):
        frame = mdanalysis_to_frame_data(empty_universe)


def test_single_atom_universe(single_atom_universe):
    """
    Tests that MD analysis converter can construct a frame
    even when only positions are supplied.
    :param single_atom_universe: An MDAnalysis universe with a single default atom in it.
    """

    frame = mdanalysis_to_frame_data(single_atom_universe)
    assert frame.particle_count == 1
    assert np.allclose(frame.particle_positions, [0, 0, 0])
    assert frame.residue_count == 1
    with pytest.raises(MissingDataError):
        _ = frame.bond_pairs


@pytest.mark.parametrize(
    "universe_attribute, mda_attribute, frame_field",
    ALL_MDA_ATTRIBUTES
)
def test_framedata_to_mda_attributes(universe_attribute, mda_attribute, frame_field, frame_data_and_universe):
    frame, universe = frame_data_and_universe
    new_universe = frame_data_to_mdanalysis(frame)
    new_universe_attribute = _get_mda_attribute(new_universe, universe_attribute, mda_attribute)
    universe_attribute = _get_mda_attribute(universe, universe_attribute, mda_attribute)
    assert all(a == b for a, b in zip(new_universe_attribute, universe_attribute))


def test_framedata_to_mda_positions(frame_data_and_universe):
    frame, universe = frame_data_and_universe
    new_universe = frame_data_to_mdanalysis(frame)
    assert np.allclose(new_universe.atoms.positions, universe.atoms.positions)


@pytest.mark.parametrize("mda_attribute, frame_field",
                         [(key, value) for key, value in MDANALYSIS_COUNTS_TO_FRAME_DATA.items()]
                         )
def test_framedata_to_mdanalysis_counts(mda_attribute, frame_field, frame_data_and_universe):
    frame, universe = frame_data_and_universe
    new_universe = frame_data_to_mdanalysis(frame)
    assert len(getattr(universe, mda_attribute)) == len(getattr(new_universe, mda_attribute))


def test_framedata_to_mda_missing_data(frame_data_and_universe):
    """
    Tests that an MDA universe can be constructed from just particle count and particle positions,
    with the resulting structure conforming to expectations.
    """
    frame, universe = frame_data_and_universe
    new_frame = FrameData()
    new_frame.particle_count = frame.particle_count
    new_frame.particle_positions = frame.particle_positions

    new_universe = frame_data_to_mdanalysis(new_frame)
    assert np.allclose(new_universe.atoms.positions, universe.atoms.positions)
    assert len(new_universe.residues) == 1
    assert len(new_universe.segments) == 1
    with pytest.raises(AttributeError):
        _ = new_universe.residues.resnames


def test_multiletter_element_symbols(metal_universe):
    frame = FrameData()
    add_mda_topology_to_frame_data(metal_universe, frame)
    assert frame.particle_elements == [11, 17, 26]
