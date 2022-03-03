# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.

"""
Methods for transmitting a simulation frame from ASE.
"""

from typing import Callable
from ase import Atoms  # type: ignore
from ase.lattice.cubic import FaceCenteredCubic
from ase.md import Langevin
from narupa.trajectory import FrameServer, FramePublisher
from narupa.ase import ase_to_frame_data


def send_ase_frame(
        ase_atoms: Atoms, frame_publisher: FramePublisher) -> Callable[[], None]:
    """
    Hook to transmit the current state of an ASE Atoms as a frame.

    :param ase_atoms: ASE :class:`Atoms`  object from which to extract frame.
    :param frame_publisher: The frame publisher on which to send the produced
        :class:`narupa.trajectory.FrameData`.

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
        frame = ase_to_frame_data(ase_atoms)
        frame_publisher.send_frame(frame_index, frame)
        frame_index += 1

    return send
