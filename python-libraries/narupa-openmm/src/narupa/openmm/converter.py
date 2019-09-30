from simtk.openmm import State
from simtk.openmm.app.topology import Topology
from simtk.unit import nanometer

from narupa.trajectory import FrameData


def add_openmm_state_to_frame_data(data: FrameData, state: State):
    positions = state.getPositions()
    box_vectors = state.getPeriodicBoxVectors()
    data.particle_positions = positions.value_in_unit(nanometer)
    data.particle_count = len(positions)
    data.box_vectors = box_vectors.value_in_unit(nanometer)


def add_openmm_topology_to_frame_data(data: FrameData, topology: Topology):
    data.residue_names = [residue.name for residue in topology.residues()]
    data.residue_ids = [residue.id for residue in topology.residues()]
    data.residue_chains = [residue.chain.index for residue in topology.residues()]
    data.residue_count = topology.getNumResidues()

    data.chain_names = [chain.name for chain in topology.chains()]
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
        bonds.append((bond[0].index, bond[1].index))

    data.particle_names = atom_names
    data.particle_elements = elements
    data.particle_residues = residue_indices
    data.particle_count = topology.getNumAtoms()

    data.bonds = bonds

    data.particle_count = len(atom_names)


def openmm_to_frame_data(*, state=None, topology=None) -> FrameData:
    data = FrameData()
    if state is not None:
        add_openmm_state_to_frame_data(data, state)
    if topology is not None:
        add_openmm_topology_to_frame_data(data, topology)
    return data
