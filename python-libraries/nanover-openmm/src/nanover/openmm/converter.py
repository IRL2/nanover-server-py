"""
Module providing conversion methods between NanoVer and OpenMM.
"""

from typing import Optional


try:
    # Mypy does not find State in the C module; we ignore the error
    from openmm import State  # type: ignore[attr-defined]
except (ImportError, ModuleNotFoundError):
    from simtk.openmm import State
try:
    from openmm.app.topology import Topology
except (ImportError, ModuleNotFoundError):
    from simtk.openmm.app.topology import Topology
from openmm.unit import kilojoule_per_mole, picosecond
from nanover.trajectory import FrameData


def add_openmm_state_to_frame_data(
    data: FrameData,
    state: State,
    include_positions: bool = True,
    include_energies: bool = True,
    end: Optional[int] = None,
) -> None:
    """
    Adds the OpenMM state information to the given :class:`FrameData`, including
    positions and periodic box vectors.

    :param data: NanoVer :class:`FrameData` to add state information to.
    :param state: OpenMM :class:`State` from which to extract state information.
    :param include_positions: If ``True``, the particle positions are read from
        the state and included in the frame.
    """
    # Here, we count of the fact that OpenMM default length unit is the
    # nanometer. By doing this assumption, we avoid arrays being copied during
    # unit conversion.
    if include_positions:
        positions = state.getPositions(asNumpy=True)[slice(None, end)]
        data.particle_positions = positions
    if include_energies:
        potential_energy = state.getPotentialEnergy().value_in_unit(kilojoule_per_mole)
        kinetic_energy = state.getKineticEnergy().value_in_unit(kilojoule_per_mole)
        total_energy = potential_energy + kinetic_energy
        data.kinetic_energy = kinetic_energy
        data.potential_energy = potential_energy
        data.total_energy = total_energy
    box_vectors = state.getPeriodicBoxVectors(asNumpy=True)
    data.box_vectors = box_vectors
    simulation_time = state.getTime().value_in_unit(picosecond)
    data.simulation_time = simulation_time


def add_openmm_topology_to_frame_data(data: FrameData, topology: Topology, end: Optional[int] = None) -> None:
    """
    Adds the OpenMM topology information to the given :class:`FrameData`,
    including residue, chain, atomic and bond information.

    :param data: :class:`FrameData` to add topology information to.
    :param topology: OpenMM :class:`Topology` from which to extract information.
    :param end: Atom index where to stop the selection of atoms to include in
        the frame. The atom with the given is excluded.
    """
    if end is None:
        end = topology.getNumAtoms()

    residues = [
        residue for residue in topology.residues()
        if max(atom.index for atom in residue.atoms()) < end
    ]
    data.residue_names = [residue.name for residue in residues]
    data.residue_ids = [residue.id for residue in residues]
    data.residue_chains = [residue.chain.index for residue in residues]
    data.residue_count = len(residues)

    chains = [
        chain for chain in topology.chains()
        if max(atom.index for atom in chain.atoms()) < end
    ]
    data.chain_names = [chain.id for chain in chains]
    data.chain_count = len(chains)

    atom_names = []
    elements = []
    residue_indices = []
    bonds = []

    for atom in (atom for atom in topology.atoms() if atom.index < end):
        atom_names.append(atom.name)
        elements.append(atom.element.atomic_number)
        residue_indices.append(atom.residue.index)

    for bond in topology.bonds():
        if bond[0].index < end and bond[1].index < end:
            bonds.append((bond[0].index, bond[1].index))

    data.particle_names = atom_names
    data.particle_elements = elements
    data.particle_residues = residue_indices
    data.particle_count = len(atom_names)

    data.bond_pairs = bonds


def openmm_to_frame_data(
    *,
    state: Optional[State] = None,
    topology: Optional[Topology] = None,
    include_positions: bool = True,
    include_energies: bool = True,
        end: Optional[int] = None
) -> FrameData:
    """
    Converts the given OpenMM state and topology objects into a Narupa :class:`FrameData`.

    Both fields are optional. For performance reasons, it is best to construct
    a Narupa :class:`FrameData` once with topology information, and from then
    on just update the state, as that will result in less data being transmitted.

    :param state: An optional OpenMM :class:`State` from which to extract
        state data.
    :param topology: An optional OpenMM :class:`Topology` from which to extract
        topological information.
    :param include_positions: If ``True``, the particle positions are read from
        the state and included in the frame.
    :return: A :class:`FrameData` with the state and topology information
        provided added to it.
    """
    data = FrameData()
    if state is not None:
        add_openmm_state_to_frame_data(data, state, include_positions, include_energies, end=end)
    if topology is not None:
        add_openmm_topology_to_frame_data(data, topology, end=end)
    return data

