from simtk.openmm.app.topology import Topology
from simtk.unit import Quantity, nanometer

from narupa.trajectory import FrameData


def add_openmm_positions_to_frame_data(data: FrameData, positions: Quantity):
    data.positions = positions.value_in_unit(nanometer)


def add_openmm_topology_to_frame_data(data: FrameData, topology: Topology):
    data.arrays['residue.id'] = [residue.name for residue in topology.residues()]
    data.arrays['residue.chain'] = [residue.chain.index for residue in topology.residues()]

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

    data.arrays['atom.id'] = atom_names
    data.elements = elements
    data.arrays['atom.residue'] = residue_indices
    data.bonds = bonds


def openmm_to_frame_data(*, positions=None, topology=None) -> FrameData:
    data = FrameData()
    if positions is not None:
        add_openmm_positions_to_frame_data(data, positions)
    if topology is not None:
        add_openmm_topology_to_frame_data(data, topology)
    return data
