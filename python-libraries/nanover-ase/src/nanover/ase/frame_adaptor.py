"""
Methods for transmitting a simulation frame from ASE.
"""

from typing import Callable
from ase import Atoms  # type: ignore

from nanover.ase.converter import (
    ase_atoms_to_topology_frame,
    ase_atoms_to_regular_frame,
)
from nanover.trajectory import FramePublisher


def send_ase_frame(
    ase_atoms: Atoms,
    frame_publisher: FramePublisher,
    include_velocities=False,
    include_forces=False,
) -> Callable[[], None]:
    """
    Hook to transmit the current state of an ASE Atoms as a frame.

    :param ase_atoms: ASE :class:`Atoms`  object from which to extract frame.
    :param frame_publisher: The frame publisher on which to send the produced
        :class:`nanover.trajectory.FrameData`.

    When attached to an ASE simulation, such as a :class:`Langevin` dynamics
    simulation, this method will be called to send the frame on the given
    :class:`FrameServer`.

    Example
    =======

    >>> frame_server = FrameServer(address="localhost", port=54321)
    >>> atoms = FaceCenteredCubic(directions=[[1, 0, 0], [0, 1, 0], [0, 0, 1]],
    ...                           symbol="Cu", size=(2, 2, 2), pbc=True)
    >>> dynamics = Langevin(atoms, timestep=0.5, temperature_K=300, friction=1.0)
    >>> dynamics.attach(send_ase_frame(atoms, frame_publisher), interval=2)
    """

    frame_index = 0

    def send():
        nonlocal frame_index
        if frame_index == 0:
            frame = ase_atoms_to_topology_frame(
                ase_atoms,
                include_velocities=include_velocities,
                include_forces=include_forces,
            )
        else:
            frame = ase_atoms_to_regular_frame(
                ase_atoms,
                include_velocities=include_velocities,
                include_forces=include_forces,
            )
        frame_publisher.send_frame(frame_index, frame)
        frame_index += 1

    return send
