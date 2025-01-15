"""
Module providing conversion methods between NanoVer and OpenMM.
"""

from typing import Optional

from openmm import State
from openmm.app.topology import Topology
from openmm.unit import kilojoule_per_mole, picosecond
from nanover.trajectory import FrameData


def add_openmm_state_to_frame_data(
    data: FrameData,
    state: State,
    include_positions=True,
    include_energies=True,
    include_velocities=False,
    include_forces=False,
    state_excludes_imd=False,
) -> None:
    """
    Adds the OpenMM state information to the given :class:`FrameData`, including
    positions and periodic box vectors.

    :param data: NanoVer :class:`FrameData` to add state information to.
    :param state: OpenMM :class:`State` from which to extract state information.
    :param include_positions: If ``True``, the particle positions are read from
        the state and included in the frame.
    :param include_energies: If ``True``, the kinetic and potential energies
        are read from the state and included in the frame.
    :param include_velocities: If ``True``, the particle velocities are read
        from the state and included in the frame.
    :param include_forces: If ``True``, the particle forces are read from the
        state and included in the frame.
    :param state_excludes_imd: Should be ``True`` if the state excludes the
        IMD force contribution.
    """
    # Here, we count of the fact that OpenMM default length unit is the
    # nanometer. By doing this assumption, we avoid arrays being copied during
    # unit conversion.
    if include_positions:
        positions = state.getPositions(asNumpy=True)
        data.particle_positions = positions
    if include_energies:
        potential_energy = state.getPotentialEnergy().value_in_unit(kilojoule_per_mole)
        kinetic_energy = state.getKineticEnergy().value_in_unit(kilojoule_per_mole)
        data.kinetic_energy = kinetic_energy
        data.potential_energy = potential_energy
    if include_velocities:
        velocities = state.getVelocities(asNumpy=True)
        data.particle_velocities = velocities
    if include_forces:
        forces = state.getForces(asNumpy=True)
        if state_excludes_imd:
            data.particle_forces_system = forces
        else:
            data.particle_forces = forces
    box_vectors = state.getPeriodicBoxVectors(asNumpy=True)
    data.box_vectors = box_vectors
    simulation_time = state.getTime().value_in_unit(picosecond)
    data.simulation_time = simulation_time


def add_openmm_topology_to_frame_data(data: FrameData, topology: Topology) -> None:
    """
    Adds the OpenMM topology information to the given :class:`FrameData`,
    including residue, chain, atomic and bond information.

    :param data: :class:`FrameData` to add topology information to.
    :param topology: OpenMM :class:`Topology` from which to extract information.
    """

    data.residue_names = [residue.name for residue in topology.residues()]
    data.residue_ids = [residue.id for residue in topology.residues()]
    data.residue_chains = [residue.chain.index for residue in topology.residues()]
    data.residue_count = topology.getNumResidues()

    data.chain_names = [chain.id for chain in topology.chains()]
    data.chain_count = topology.getNumChains()

    atom_names = []
    elements = []
    residue_indices = []
    bonds = []

    for atom in topology.atoms():
        atom_names.append(atom.name)
        elements.append(atom.element.atomic_number)
        residue_indices.append(atom.residue.index)

    for bond in topology.bonds():
        bonds.append([bond[0].index, bond[1].index])

    data.particle_names = atom_names
    data.particle_elements = elements
    data.particle_residues = residue_indices
    data.particle_count = topology.getNumAtoms()
    data.bond_pairs = bonds


def openmm_to_frame_data(
    *,
    state: Optional[State] = None,
    topology: Optional[Topology] = None,
    include_positions=True,
    include_energies=True,
    include_velocities=False,
    include_forces=False,
    state_excludes_imd=False,
) -> FrameData:
    """
    Converts the given OpenMM state and topology objects into a NanoVer :class:`FrameData`.

    Both fields are optional. For performance reasons, it is best to construct
    a NanoVer :class:`FrameData` once with topology information, and from then
    on just update the state, as that will result in less data being transmitted.

    :param state: An optional OpenMM :class:`State` from which to extract
        state data.
    :param topology: An optional OpenMM :class:`Topology` from which to extract
        topological information.
    :param include_positions: If ``True``, the particle positions are read from
        the state and included in the frame.
    :param include_energies: If ``True``, the kinetic and potential energies
        are read from the state and included in the frame.
    :param include_velocities: If ``True``, the particle velocities are read
        from the state and included in the frame.
    :param include_forces: If ``True``, the particle forces are read from the
        state and included in the frame.
    :return: A :class:`FrameData` with the state and topology information
        provided added to it.
    """
    data = FrameData()
    if state is not None:
        add_openmm_state_to_frame_data(
            data,
            state,
            include_positions,
            include_energies,
            include_velocities,
            include_forces,
            state_excludes_imd,
        )
    if topology is not None:
        add_openmm_topology_to_frame_data(data, topology)
    return data
