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
        # generate topology frame using OpenMM converter.
        if frame_index == 0:
            frame = openmm_ase_atoms_to_topology_frame(ase_atoms, **kwargs)
        # from then on, just send positions and state.
        else:
            frame = openmm_ase_atoms_to_regular_frame(ase_atoms, **kwargs)
        frame_publisher.send_frame(frame_index, frame)
        frame_index += 1

    return send


def openmm_ase_atoms_to_regular_frame(
    ase_atoms: Atoms,
    **kwargs,
):
    frame = ase_to_frame_data(
        ase_atoms,
        topology=False,
        **kwargs,
    )
    return frame


def openmm_ase_atoms_to_topology_frame(
    ase_atoms: Atoms,
    **kwargs,
):
    imd_calculator = ase_atoms.calc
    topology = imd_calculator.calculator.topology
    frame = openmm_ase_atoms_to_regular_frame(ase_atoms, **kwargs)
    add_openmm_topology_to_frame_data(frame, topology)
    return frame
