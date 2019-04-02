"""
Convert data structures from OpenMM to Narupa.
"""
# DRPC Data classes have members that are not found by pylint. We ignore the
# pylint error to avoid noisy reports.
# pylint: disable=no-member

from simtk.openmm.app.topology import Topology
from simtk.unit.quantity import Quantity

# The prefixed units are programmatically added to the simtk.unit module, thus
# there are not found by pyint or PyCharm.
from simtk.unit import nanometer  # pylint: disable=no-name-in-module

from narupa.protocol.trajectory.frame_pb2 import FrameData
from narupa.protocol.topology.topology_pb2 import TopologyData


def openmm_positions_to_frame_data(positions: Quantity) -> FrameData:
    """
    Convert OpenMM positions to GRPC ready Narupa positions.

    :param positions: Positions in OpenMM format; for instance as returned by
        :meth:`simtk.openmm.State.getPositions`.
    :return: The positions in a format ready to be sent to GRPC.
    """
    data = FrameData()

    array = data.arrays['atom.position'].float_values.values

    floats = [
        value
        for position in positions
        for value in position.value_in_unit(nanometer)
    ]
    array.extend(floats)

    return data


def openmm_topology_to_topology_data(topology: Topology) -> TopologyData:
    """
    Convert an OpenMM topology to a GRPC ready Narupa topology.

    :param topology: An instance of OpenMM :class:`Topology`.
    :return: The topology in a format ready to be sent to GRPC.
    """
    data = TopologyData()

    data.arrays['residue.id'].string_values.values.extend(
        [residue.name for residue in topology.residues()]
    )
    data.arrays['residue.chain'].index_values.values.extend(
        [residue.chain.index for residue in topology.residues()]
    )

    for atom in topology.atoms():
        data.arrays['atom.id'].string_values.values.append(atom.name)
        data.arrays['atom.element'].index_values.values.append(atom.element.atomic_number)
        data.arrays['atom.residue'].index_values.values.append(atom.residue.index)

    for bond in topology.bonds():
        data.arrays['bond'].index_values.values.append(bond[0].index)
        data.arrays['bond'].index_values.values.append(bond[1].index)

    return data
