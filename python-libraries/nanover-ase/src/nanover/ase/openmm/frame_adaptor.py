from ase import Atoms

from nanover.ase import ase_to_frame_data
from nanover.openmm.converter import add_openmm_topology_to_frame_data
from nanover.trajectory import FramePublisher


def openmm_ase_frame_adaptor(
    ase_atoms: Atoms,
    frame_publisher: FramePublisher,
    **kwargs,
):
    """
    Generates and sends frames for a simulation using an :class: OpenMMCalculator.
    """

    frame_index = 0

    def send():
        nonlocal frame_index
        include_topology = frame_index == 0
        frame_data = openmm_ase_atoms_to_frame_data(
            ase_atoms, topology=include_topology, **kwargs
        )
        frame_publisher.send_frame(frame_index, frame_data)
        frame_index += 1

    return send


def openmm_ase_atoms_to_frame_data(
    ase_atoms: Atoms,
    *,
    topology=False,
    **kwargs,
):
    frame_data = ase_to_frame_data(
        ase_atoms,
        topology=False,
        **kwargs,
    )

    if topology:
        imd_calculator = ase_atoms.calc
        topology = imd_calculator.calculator.topology
        add_openmm_topology_to_frame_data(frame_data, topology)

    return frame_data
