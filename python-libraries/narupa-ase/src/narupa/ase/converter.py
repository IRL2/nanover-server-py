from typing import Iterable

from ase import Atoms
import itertools
import numpy as np

from narupa.protocol.trajectory import FrameData


ANG_TO_NM = 0.1
NM_TO_ANG = 1.0 / ANG_TO_NM

# from https://chem.libretexts.org/Bookshelves/Physical_and_Theoretical_Chemistry_Textbook_Maps/Supplemental_Modules_(Physical_and_Theoretical_Chemistry)/Chemical_Bonding/Fundamentals_of_Chemical_Bonding/Covalent_Bond_Distance%2C_Radius_and_van_der_Waals_Radius
# TODO make a helper class with complete coverage of this stuff.
ATOM_RADIUS_ANG = {
    "H": 1.2,
    "B": 2.08,
    "C": 1.85,
    "N": 1.54,
    "O": 1.40,
    "F": 1.35,
    "Cl": 1.80,
    "Br": 1.95,
    "I": 2.15,
    "He": 0.99,
    "Cu": 2.38,
    "Pt": 2.32}

def ase_to_framedata(ase_atoms: Atoms, positions=True, topology=True, state=True) -> FrameData:
    data = FrameData()

    if topology:
        data.arrays['residue.id'].string_values.values.extend(["ASE"])
        data.arrays['residue.chain'].index_values.values.extend([0])

        for index, atom in enumerate(ase_atoms):
            data.arrays['atom.id'].string_values.values.append(str(atom.index))
            data.arrays['atom.element'].index_values.values.append(atom.number)
            data.arrays['atom.residue'].index_values.values.append(0)

        bonds = generate_bonds(ase_atoms)
        for bond in bonds:
            data.arrays['bond'].index_values.values.append(bond[0])
            data.arrays['bond'].index_values.values.append(bond[1])

    if positions:
        array = data.arrays['atom.position'].float_values.values
        floats = ase_atoms.get_positions().flatten() * ANG_TO_NM
        array[:] = floats

    if state:
        data.values['energy.potential'].number_value = ase_atoms.get_potential_energy()
        data.values['energy.kinetic'].number_value = ase_atoms.get_kinetic_energy()

    return data


def get_radius_of_element(symbol):
    """
    Gets the radius of an atom (Angstrom)
    :param symbol:
    :return:
    """
    return ATOM_RADIUS_ANG[symbol]


def bond_threshold(radii: Iterable):
    """
    Returns the distance threshold for a bond between a given pair of radii.
    :param radii: Pair or radii.
    :return: Distance threshold indicating the distance at which a bond is formed.
    """
    return 0.6 * sum(radii)


def generate_bonds(atoms: Atoms):
    bonds = []
    for pair in itertools.combinations(atoms, 2):
        distance = np.linalg.norm(pair[0].position - pair[1].position)
        radii = [get_radius_of_element(atom.symbol) for atom in pair]
        if distance < bond_threshold(radii):
            bonds.append([atom.index for atom in pair])
    return bonds
