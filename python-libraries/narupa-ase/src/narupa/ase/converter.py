from typing import Iterable, Optional

from ase import Atoms, Atom
import itertools
import numpy as np

from narupa.trajectory import FrameData


ANG_TO_NM = 0.1
NM_TO_ANG = 1.0 / ANG_TO_NM
KJMOL_TO_EV = 0.01036427
EV_TO_KJMOL = 1.0 / KJMOL_TO_EV

# from https://chem.libretexts.org/Bookshelves/Physical_and_Theoretical_Chemistry_Textbook_Maps/Supplemental_Modules_(Physical_and_Theoretical_Chemistry)/Chemical_Bonding/Fundamentals_of_Chemical_Bonding/Covalent_Bond_Distance%2C_Radius_and_van_der_Waals_Radius
# TODO make a helper class with complete coverage of this stuff.
ATOM_RADIUS_ANG = {
    "H": 1.2,
    "He": 1.4,
    "Li": 1.82,
    "Be": 1.53,
    "B": 1.92,
    "C": 1.70,
    "N": 1.55,
    "O": 1.52,
    "F": 1.47,
    "Ne": 1.54,
    "Na": 2.27,
    "Mg": 1.73,
    "Al": 1.84,
    "Si": 2.1,
    "P": 1.8,
    "S": 1.8,
    "Cl": 1.75,
    "Ar": 1.88,
    "K": 2.75,
    "Ca": 2.31,
    "Ni": 1.63,
    "Cu": 1.40,
    "Zn": 1.39,
    "Ga": 1.87,
    "Ge": 2.11,
    "As": 1.85,
}

def ase_to_frame_data(ase_atoms: Atoms, positions=True, topology=True, state=True) -> FrameData:
    """
    Constructs a Narupa frame from the state if an ASE simulation.
    :param ase_atoms: The ASE atoms object representing the state of the simulation to send.
    :param positions: Whether to add positions to the frame.
    :param topology: Whether to generate the current state of the topology and add it to the frame.
    :param state: Whether to add additional state information such as energies.
    :return: Narupa frame.
    """
    data = FrameData()
    if positions:
        add_ase_positions_to_frame_data(data, ase_atoms.get_positions())
    if topology:
        add_ase_topology_to_frame_data(data, ase_atoms)
    if state:
        add_ase_state_to_frame_data(data, ase_atoms)
    return data




def frame_data_to_ase(frame_data: FrameData, positions=True, topology=True, ase_atoms=None) -> Atoms:
    if ase_atoms is None:
        ase_atoms = Atoms()
    if topology:
        ase_atoms = Atoms()
        add_frame_data_topology_to_ase(frame_data, ase_atoms)
    if positions:
        add_frame_data_positions_to_ase(frame_data, ase_atoms)
    return ase_atoms

def add_frame_data_topology_to_ase(data: FrameData, atoms: Atoms):
    for element in data.particle_elements:
        atoms.append(Atom(symbol=element))

def add_frame_data_positions_to_ase(frame_data, ase_atoms):
    ase_atoms.set_positions(np.array(frame_data.particle_positions) * NM_TO_ANG)

def add_ase_positions_to_frame_data(data: FrameData, positions: np.array ):
    """
    Adds ASE positions to the frame data, converting to nanometers.
    :param data:
    :param positions:
    """
    data.particle_positions = positions * ANG_TO_NM

def add_ase_topology_to_frame_data(data: FrameData, ase_atoms: Atoms):
    """
    Generates a topology for the current state of the atoms and adds it to the frame.
    :param data:
    :param ase_atoms:
    """
    data.arrays['residue.id'] = ["ASE"]
    data.arrays['residue.chain'] = [0]

    atom_names = []
    elements = []
    residue_ids = []
    for index, atom in enumerate(ase_atoms):
        atom_names.append(str(atom.index))
        elements.append(atom.number)
        residue_ids.append(0)

    data.arrays['atom.id'] = atom_names
    data.particle_elements = elements
    data.arrays['atom.residue'] = residue_ids

    bonds = generate_bonds(ase_atoms)
    data.bonds = bonds

def add_ase_state_to_frame_data(data: FrameData, ase_atoms: Atoms):
    """
    Adds simulaton state information to the frame.
    :param data:
    :param ase_atoms:
    :return:
    """
    # get the energy directly, without performing a recalculation.
    energy = ase_atoms.get_calculator().get_property('energy', allow_calculation=False)
    if energy is not None:
        data.values['energy.potential'] = energy
    data.values['energy.kinetic'] = ase_atoms.get_kinetic_energy()

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
