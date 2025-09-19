"""
Module containing methods for converting between ASE simulations consisting of
:class:`Atoms`  and the NanoVer :class:`FrameData` object for transmission to
NanoVer clients.
"""

from typing import Iterable

from ase import Atoms, Atom  # type: ignore
from ase.units import fs as fs_in_ase_time_unit
import itertools
import numpy as np
import numpy.typing as npt

from nanover.ase.imd_calculator import ImdCalculator
from nanover.trajectory import FrameData

ANG_TO_NM = 0.1
NM_TO_ANG = 1.0 / ANG_TO_NM
KJMOL_TO_EV = 0.01036427
EV_TO_KJMOL = 1.0 / KJMOL_TO_EV
PS_TO_ASE_TIME_UNIT = fs_in_ase_time_unit * 1e3
ASE_TIME_UNIT_TO_PS = 1.0 / PS_TO_ASE_TIME_UNIT

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


def ase_atoms_to_frame_data(
    ase_atoms: Atoms,
    *,
    topology: bool,
    **kwargs,
) -> FrameData:
    return ase_to_frame_data(ase_atoms, topology=topology, **kwargs)


def ase_to_frame_data(
    ase_atoms: Atoms,
    positions=True,
    topology=True,
    state=True,
    box_vectors=True,
    generate_bonds=True,
    include_velocities=False,
    include_forces=False,
) -> FrameData:
    """
    Constructs a NanoVer frame from the state of the atoms in an ASE simulation.

    :param ase_atoms: The ASE atoms object representing the state of the
        simulation to send.
    :param positions: Whether to add positions to the frame.
    :param topology: Whether to generate the current state of the topology and
        add it to the frame.
    :param state: Whether to add additional state information such as energies.
    :param box_vectors: Whether to add the box vectors to the frame data.
    :param generate_bonds: Whether to generate bonds for the topology.
    :param include_velocities: Whether to includes per particle velocities.
    :param include_forces: Whether to include per particle forces.
    :return: NanoVer frame.

    :raises: AttributeError Raised if state is `True`, and `ase_atoms` has no
        calculator attached.

    Example
    =======

    >>> atoms = Atoms('CO', positions=[(0, 0, 0), (0, 0, 1.1)], cell=[2, 2, 2])
    >>> frame = ase_to_frame_data(atoms, state=False)
    >>> frame.particle_count
    2
    >>> frame.bonds
    [[0, 1]]
    >>> frame.particle_elements
    [6, 8]

    """
    data = FrameData()
    if positions:
        add_ase_positions_to_frame_data(data, ase_atoms.get_positions(wrap=False))
    if topology:
        add_ase_topology_to_frame_data(data, ase_atoms, generate_bonds=generate_bonds)
    if state:
        add_ase_state_to_frame_data(data, ase_atoms)
    if box_vectors:
        add_ase_box_vectors_to_frame_data(data, ase_atoms)
    if include_velocities:
        add_ase_velocities_to_frame_data(data, ase_atoms)
    if include_forces:
        add_ase_forces_to_frame_data(data, ase_atoms)

    return data


def frame_data_to_ase(
    frame_data: FrameData,
    positions: bool = True,
    topology: bool = True,
    ase_atoms: Atoms | None = None,
) -> Atoms:
    """
    Constructs an ASE :class:`Atoms` object from a NanoVer :class:`FrameData`.

    :param frame_data: The NanoVer :class:`FrameData`.
    :param positions: Whether to add positions to the ASE atoms.
    :param topology: Whether to add topology information within the frame data
        to ASE.
    :param ase_atoms: Optional ASE :class:`Atoms` object, which will have its
        positions replaced. If the flag `topology` is set, then a new object
        will still be constructed.
    :return: ASE Atoms updated or created with the data contained in the
        NanoVer frame.

    Example:
    ========

    >>> frame = FrameData()
    >>> frame.particle_elements = [6, 8]
    >>> frame.particle_positions = [[0,0,0], [0,0, 0.11]]
    >>> atoms = frame_data_to_ase(frame)
    >>> atoms.symbols
    Symbols('CO')

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
    Adds frame data topology information to ASE :class:`Atoms`.

    Since ASE :class:`Atoms` do not have a concept of bonds, this just adds
    particle elements.

    :param frame_data: :class:`FrameData` from which to extract topology.
    :param atoms: ASE :class:`Atoms` to add element data to.
    """
    for element in frame_data.particle_elements:
        atoms.append(Atom(symbol=element))
    frame_data.particle_count = len(atoms)


def add_frame_data_positions_to_ase(frame_data, ase_atoms):
    """
    Adds frame data particle positions to ASE atoms, converting to angstroms.

    :param frame_data: :class:`FrameData` from which to extract positions.
    :param ase_atoms: ASE :class:`Atoms` to add particle positions to.
    """
    ase_atoms.set_positions(np.array(frame_data.particle_positions) * NM_TO_ANG)


def add_ase_positions_to_frame_data(data: FrameData, positions: npt.NDArray):
    """
    Adds ASE positions to the frame data, converting to nanometers.

    :param data: :class:`FrameData` to add atom positions to.
    :param positions: Array of atomic positions, in angstroms.
    """
    data.particle_positions = positions * ANG_TO_NM


def add_ase_velocities_to_frame_data(data: FrameData, ase_atoms: Atoms):
    """
    Adds ASE velocities to the frame data, converting to nanometers per picosecond.

    :param data: :class:`FrameData` to add atom velocities to.
    :param ase_atoms: ASE :class:`Atoms` to add particle positions to.
    """
    data.particle_velocities = ase_atoms.get_velocities() * (
        ANG_TO_NM / ASE_TIME_UNIT_TO_PS
    )


def add_ase_forces_to_frame_data(data: FrameData, ase_atoms: Atoms):
    """
    Adds ASE forces to the frame data, converting to kJ mol-1 per nanometer. If the ASE
    calculator is an ImdCalculator, removes the iMD forces from the ASE forces to deliver
    the system forces (iMD forces delivered separately elsewhere).

    :param data: :class:`FrameData` to add atom forces to.
    :param ase_atoms: ASE :class:`Atoms` to add particle positions to.
    """

    data.particle_forces_system = ase_atoms.get_forces() * (EV_TO_KJMOL / ANG_TO_NM)
    if isinstance(ase_atoms.calc, ImdCalculator):
        data.particle_forces_system -= ase_atoms.calc.results["interactive_forces"] * (
            EV_TO_KJMOL / ANG_TO_NM
        )


def add_ase_box_vectors_to_frame_data(data: FrameData, ase_atoms: Atoms):
    """
    Adds the periodic box vectors from the given ASE :class:`Atoms`
    object to the given :class:`FrameData`.

    :param data: :class:`FrameData` upon which to add periodic box vectors.
    :param ase_atoms: :class:`Atoms` from which to extract periodic box vectors.
    """
    box_vectors = ase_atoms.cell.copy() * ANG_TO_NM
    data.box_vectors = box_vectors


def add_ase_topology_to_frame_data(
    frame_data: FrameData, ase_atoms: Atoms, generate_bonds=True
):
    """
    Generates a topology for the current state of the atoms and adds it to the frame.

    Since ASE atoms have no concept of bonds, they are generated using distance
    criteria.

    :param frame_data: Frame data to add topology information to.
    :param ase_atoms: ASE atoms to extract topology information from.
    """
    # TODO it would be nice to do dynamic molecule/chain detection here.
    frame_data.residue_names = ["ASE"]
    frame_data.residue_chains = [0]
    frame_data.residue_count = 1
    frame_data.residue_ids = ["1"]

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

    if generate_bonds:
        bonds = generate_bonds_from_ase(ase_atoms)
        frame_data.bond_pairs = bonds


def add_ase_state_to_frame_data(frame_data: FrameData, ase_atoms: Atoms):
    """
    Adds simulation state information to the frame,
    consisting of the potential energy and kinetic energy of the
    molecular system. If the ASE calculator is an ImdCalculator,
    removes the iMD energy from the ASE potential energy to deliver
    the system potential energy (iMD energy delivered separately
    elsewhere).

    :param frame_data: Frame data to add ASE state information to.
    :param ase_atoms: The ASE atoms from which to extract state information.
    """
    # get the energy directly, without performing a recalculation.
    try:
        energy = ase_atoms.calc.get_property("energy", allow_calculation=False)
    except AttributeError:
        raise AttributeError("No calculator in atoms, so energy cannot be computed")
    if energy is not None:
        frame_data.potential_energy = energy * EV_TO_KJMOL
        # Subtract iMD energy from total potential energy to obtain system potential energy
        if isinstance(ase_atoms.calc, ImdCalculator):
            frame_data.potential_energy -= (
                ase_atoms.calc.get_property(
                    "interactive_energy", allow_calculation=False
                )
                * EV_TO_KJMOL
            )
    frame_data.kinetic_energy = ase_atoms.get_kinetic_energy() * EV_TO_KJMOL


def get_radius_of_element(symbol: str, default=1.0):
    """
    Gets the radius of an atom in Angstroms.

    :param symbol: The chemical symbol representing the element.
    :param default: Default radius to use if the radius for the given chemical
        symbol is not known.
    :return: The VDW radius of the atom in angstroms.
    """
    return ATOM_RADIUS_ANG.get(symbol, default)


def _bond_threshold(radii: Iterable):
    """
    Returns the distance threshold, in Angstroms, for a bond between a given pair of radii.

    :param radii: Pair or radii.
    :return: Distance threshold indicating the distance at which a bond is formed.
    """
    return 0.6 * sum(radii)


def generate_bonds_from_ase(atoms: Atoms):
    """
    Generates bonds for the given configuration of ASE atoms using a distance criterion.

    A bond is placed between two atoms if the distance between them is less
    than 0.6 times the VDW radii of the atoms.

    :param atoms: ASE atoms to generate bonds for.
    :return: A list of pairs of atom indexes representing bonds.

    Example
    =======

    The following example produces a bond between the two atoms in a carbon
    monoxide molecule:

    >>> co = Atoms('CO', positions=[(0, 0, 0), (0, 0, 1.1)], cell=[2, 2, 2])
    >>> generate_bonds(co)
    [[0, 1]]

    """
    bonds = []
    for pair in itertools.combinations(atoms, 2):
        distance = np.linalg.norm(pair[0].position - pair[1].position)
        radii = [get_radius_of_element(atom.symbol) for atom in pair]
        if distance < _bond_threshold(radii):
            bonds.append([atom.index for atom in pair])
    return bonds
