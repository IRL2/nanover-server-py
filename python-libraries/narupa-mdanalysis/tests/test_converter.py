import numpy as np
import pytest
from MDAnalysis import Universe
from narupa.mdanalysis.converter import mdanalysis_to_frame_data, INDEX_ELEMENT, FRAME_DATA_TO_MDANALYSIS_COUNTS, \
    ALL_MDA_ATTRIBUTES
from narupa.trajectory.frame_data import (PARTICLE_ELEMENTS, MissingDataError)

TEST_SYSTEM = "2efv_fragment.pdb"


@pytest.fixture
def universe():
    return Universe(TEST_SYSTEM, guess_bonds=True)


@pytest.fixture()
def empty_universe():
    return Universe.empty(1, trajectory=True)

@pytest.fixture()
def empty_universe_no_positions():
    return Universe.empty(1, trajectory=False)

@pytest.fixture()
def single_atom_universe(empty_universe: Universe):
    empty_universe.atoms.positions = [[0, 0, 0]]
    return empty_universe


@pytest.fixture
def frame_data(universe):
    return mdanalysis_to_frame_data(universe), universe


def test_mdanalysis_to_frame_data(universe):
    frame_data = mdanalysis_to_frame_data(universe)
    assert frame_data is not None


@pytest.mark.parametrize("universe_attribute, mda_attribute, frame_field",
                         ALL_MDA_ATTRIBUTES
                         )
def test_mdanalysis_particle_field(universe_attribute, mda_attribute, frame_field, frame_data):
    frame, universe = frame_data
    # fetches the atoms, residues or chains object, then the attribute.
    attribute = getattr(getattr(universe, universe_attribute), mda_attribute)
    field = frame.arrays[frame_field]
    if frame_field == PARTICLE_ELEMENTS:
        field = [INDEX_ELEMENT[x] for x in field]
    assert all(a == b for a, b in zip(attribute, field))


def test_mdanalysis_positions(frame_data):
    frame, universe = frame_data
    assert np.allclose(np.array(frame.particle_positions) * 10, universe.atoms.positions)


@pytest.mark.parametrize("mda_attribute, frame_field",
                         [(key, value) for key, value in FRAME_DATA_TO_MDANALYSIS_COUNTS.items()]
                         )
def test_mdanalysis_counts(mda_attribute, frame_field, frame_data):
    frame, universe = frame_data
    assert len(getattr(universe, mda_attribute)) == frame.values[frame_field]


def test_mdanalysis_bonds(frame_data):
    frame, universe = frame_data
    frame_bonds = np.array(frame.bonds)
    universe_bonds = np.array(universe.bonds.indices)
    assert np.allclose(frame_bonds, universe_bonds)


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
        _ = frame.bonds
