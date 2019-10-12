from ase.calculators.lj import LennardJones
from narupa.ase import ase_to_frame_data
import numpy as np
import pytest
from narupa.ase.converter import EV_TO_KJMOL, frame_data_to_ase
from narupa.trajectory.frame_data import MissingDataError

from .test_imd_calculator import atoms


@pytest.fixture
def atoms_lj_calc(atoms):
    calc = LennardJones()
    atoms.set_calculator(calc)
    return atoms


@pytest.fixture
def frame(atoms):
    return ase_to_frame_data(atoms, state=False)


def test_convert_atom_positions(frame):
    assert np.allclose(frame.particle_positions, [(0, 0, 0), (0, 0, 0.11)])


def test_convert_bonds(frame):
    assert np.allclose(frame.bonds, [(0, 1)])


def test_convert_atom_residues(frame):
    assert frame.particle_residues == [0, 0]


def test_convert_atom_residue_count(frame):
    assert frame.residue_count == 1


def test_convert_residue_ids(frame):
    assert frame.residue_ids == [0]


def test_convert_residue_names(frame):
    assert frame.residue_names == ['ASE']


def test_convert_residue_chains(frame):
    assert frame.residue_chains == [0]


def test_convert_chain_names(frame):
    assert frame.chain_names == ['A']


def test_convert_chain_count(frame):
    assert frame.chain_count == 1


def test_frame_positions_only(atoms):
    frame = ase_to_frame_data(atoms, positions=True, topology=False, state=False)
    with pytest.raises(MissingDataError):
        _ = frame.bonds
    with pytest.raises(MissingDataError):
        _ = frame.particle_residues
    with pytest.raises(MissingDataError):
        _ = frame.residue_chains
    with pytest.raises(MissingDataError):
        _ = frame.chain_names


def test_frame_no_positions(atoms):
    frame = ase_to_frame_data(atoms, positions=False, state=False)
    with pytest.raises(MissingDataError):
        _ = frame.particle_positions


def test_frame_state_no_calculator(atoms):
    with pytest.raises(AttributeError):
        _ = ase_to_frame_data(atoms)


def test_frame_potential_energy(atoms_lj_calc):
    energy = atoms_lj_calc.get_potential_energy()
    frame = ase_to_frame_data(atoms_lj_calc)
    assert frame.potential_energy == energy * EV_TO_KJMOL


def test_frame_potential_energy_no_calculation(atoms_lj_calc):
    """
    Tests that if no calculation has been peformed, no potential energy is produced in the frame.
    :param atoms_lj_calc: ASE atoms object with a lennard jones calculator.
    """
    frame = ase_to_frame_data(atoms_lj_calc)
    with pytest.raises(MissingDataError):
        _ = frame.potential_energy


def test_frame_kinetic_energy(atoms_lj_calc):
    atoms_lj_calc.set_velocities([(1, 1, 1), (2, 2, 2)])
    ke = atoms_lj_calc.get_kinetic_energy()
    frame = ase_to_frame_data(atoms_lj_calc)
    assert frame.kinetic_energy == ke * EV_TO_KJMOL


def test_frame_to_ase_positions(frame):
    atoms = frame_data_to_ase(frame)
    assert np.allclose(atoms.positions, [(0, 0, 0), (0, 0, 1.1)])

def test_frame_to_ase_elements(frame):
    atoms = frame_data_to_ase(frame)
    assert atoms.get_chemical_symbols() == ['C', 'O']
