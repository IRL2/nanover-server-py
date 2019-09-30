"""Narupa ASE conversion tools.

This module methods for converting between ASE simulations and Narupa frames.

"""
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


def ase_to_frame_data(
        ase_atoms: Atoms,
        positions: bool = True,
        topology: bool = True,
        state: bool = True,
        box_vectors: bool = True,
) -> FrameData:
    """
    Constructs a Narupa frame from the state of the atoms in an ASE simulation.

    :param ase_atoms: The ASE atoms object representing the state of the simulation to send.
    :param positions: Whether to add positions to the frame.
    :param topology: Whether to generate the current state of the topology and add it to the frame.
    :param state: Whether to add additional state information such as energies.
    :param box_vectors: Whether to add the box vectors to the frame data.
    :return: Narupa frame.
    """
    data = FrameData()
    if positions:
        add_ase_positions_to_frame_data(data, ase_atoms.get_positions(wrap=False))
    if topology:
        add_ase_topology_to_frame_data(data, ase_atoms)
    if state:
        add_ase_state_to_frame_data(data, ase_atoms)
    if box_vectors:
        add_ase_box_vectors_to_frame_data(data, ase_atoms)
    return data


def frame_data_to_ase(frame_data: FrameData, positions: bool = True, topology: bool = True,
                      ase_atoms: Optional[Atoms] = None) -> Atoms:
    """
    Constructs an ASE atoms object from a Narupa frame.

    :param frame_data: The Narupa frame.
    :param positions: Whether to add positions to the ASE atoms.
    :param topology: Whether to add topology information within the frame data to ASE.
    :param ase_atoms:
    :return:
    """
    if ase_atoms is None:
        ase_atoms = Atoms()
    if topology:
        ase_atoms = Atoms()
        add_frame_data_topology_to_ase(frame_data, ase_atoms)
    if positions:
        add_frame_data_positions_to_ase(frame_data, ase_atoms)
    return ase_atoms


def add_frame_data_topology_to_ase(frame_data: FrameData, atoms: Atoms):
    """
    Adds frame data topology information to ASE atoms.
    Since ASE atoms do not have a concept of bonds, this just adds
    particle elements.

    :param frame_data: Frame data from which to extract topology.
    :param atoms: ASE atoms to add element data to.
    """
    for element in frame_data.particle_elements:
        atoms.append(Atom(symbol=element))
    frame_data.particle_count = len(atoms)


def add_frame_data_positions_to_ase(frame_data, ase_atoms):
    """
    Adds frame data particle positions to ASE atoms, converting to angstroms.

    :param frame_data: Frame data from which to extract positions.
    :param ase_atoms: ASE atoms to add particle positions to.
    """
    ase_atoms.set_positions(np.array(frame_data.particle_positions) * NM_TO_ANG)


def add_ase_positions_to_frame_data(data: FrameData, positions: np.array):
    """
    Adds ASE positions to the frame data, converting to nanometers.

    :param data:
    :param positions:
    """
    data.particle_positions = positions * ANG_TO_NM
    data.particle_count = len(positions)


def add_ase_box_vectors_to_frame_data(data: FrameData, ase_atoms: Atoms):
    """
    Adds the periodic box vectors to the frame.
    """
    box_vectors = ase_atoms.cell.copy() * ANG_TO_NM
    data.box_vectors = box_vectors


def add_ase_topology_to_frame_data(frame_data: FrameData, ase_atoms: Atoms):
    """
    Generates a topology for the current state of the atoms and adds it to the frame.
    Since ASE atoms have no concept of bonds, they are generated using distance criteria.

    :param frame_data: Frame data to add topology information to.
    :param ase_atoms: ASE atoms to extract topology information from.
    """
    frame_data.residue_names = ["ASE"]
    frame_data.residue_ids = ["1"]
    frame_data.residue_chains = [0]
    frame_data.residue_count = 1

    frame_data.chain_names = ["A"]
    frame_data.chain_count = 1

    atom_names = []
    elements = []
    residue_ids = []
    for index, atom in enumerate(ase_atoms):
        atom_names.append(str(atom.index))
        elements.append(atom.number)
        residue_ids.append(0)

    frame_data.particle_names = atom_names
    frame_data.particle_elements = elements
    frame_data.particle_residues = residue_ids
    frame_data.particle_count = len(ase_atoms)

    bonds = generate_bonds(ase_atoms)
    frame_data.bonds = bonds


def add_ase_state_to_frame_data(frame_data: FrameData, ase_atoms: Atoms):
    """
    Adds simulation state information to the frame,
    consisting of the potential energy and kinetic energy.

    :param frame_data: Frame data to add ASE state information to.
    :param ase_atoms: The ASE atoms from which to extract state information.
    """
    # get the energy directly, without performing a recalculation.
    energy = ase_atoms.get_calculator().get_property('energy', allow_calculation=False)
    if energy is not None:
        frame_data.potential_energy = energy
    frame_data.kinetic_energy = ase_atoms.get_kinetic_energy()


def get_radius_of_element(symbol: str):
    """
    Gets the radius of an atom in angstroms.

    :param symbol: The chemical symbol representing the element.
    :return: The VDW radius of the atom in angstroms.
    """
    return ATOM_RADIUS_ANG[symbol]


def _bond_threshold(radii: Iterable):
    """
    Returns the distance threshold for a bond between a given pair of radii.

    :param radii: Pair or radii.
    :return: Distance threshold indicating the distance at which a bond is formed.
    """
    return 0.6 * sum(radii)


def generate_bonds(atoms: Atoms):
    """
    Generates bonds for the given configuration of ASE atoms using a distance criterion.

    A bond is placed between two atoms if the distance between them is less than 0.6 times
    the VDW radii of the atoms.

    :param atoms: ASE atoms to generate bonds for.
    :return: A list of pairs of atom indexes representing bonds.
    """
    bonds = []
    for pair in itertools.combinations(atoms, 2):
        distance = np.linalg.norm(pair[0].position - pair[1].position)
        radii = [get_radius_of_element(atom.symbol) for atom in pair]
        if distance < _bond_threshold(radii):
            bonds.append([atom.index for atom in pair])
    return bonds
