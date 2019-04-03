from simtk.openmm.app.topology import Topology

from narupa.protocol.trajectory.frame_pb2 import FrameData


def add_openmm_positions_to_frame_data(data: FrameData, positions):
    array = data.arrays['atom.position'].float_values.values

    floats = [value for position in positions for value in position._value]
    array.extend(floats)


def add_openmm_topology_to_frame_data(data: FrameData, topology: Topology):
    data.arrays['residue.id'].string_values.values.extend([residue.name for residue in topology.residues()])
    data.arrays['residue.chain'].index_values.values.extend([residue.chain.index for residue in topology.residues()])

    for atom in topology.atoms():
        data.arrays['atom.id'].string_values.values.append(atom.name)
        data.arrays['atom.element'].index_values.values.append(atom.element.atomic_number)
        data.arrays['atom.residue'].index_values.values.append(atom.residue.index)

    for bond in topology.bonds():
        data.arrays['bond'].index_values.values.append(bond[0].index)
        data.arrays['bond'].index_values.values.append(bond[1].index)


def openmm_to_frame_data(*, positions=None, topology=None) -> FrameData:
    data = FrameData()
    if positions is not None:
        add_openmm_positions_to_frame_data(data, positions)
    if topology is not None:
        add_openmm_topology_to_frame_data(data, topology)
    return data
