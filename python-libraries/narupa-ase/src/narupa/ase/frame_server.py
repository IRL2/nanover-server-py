# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.

"""
Methods for transmitting a simulation frame from ASE.
"""

from ase import Atoms
from ase.lattice.cubic import FaceCenteredCubic
from ase.md import Langevin
from narupa.trajectory import FrameServer
from narupa.ase import ase_to_frame_data


def send_ase_frame(ase_atoms: Atoms, frame_server: FrameServer):
    """
    Hook to transmit the current state of an ASE Atoms as a frame.
    :param ase_atoms: ASE :class:`Atoms`  object from which to extract frame.
    :param frame_server: The frame server on which to send the produced ~:class: narupa.trajectory.FrameData.

    When attached to an ASE simulation, such as a :class:`Langevin` dynamics simulation, this method will be
    called to send the frame on the given :class:`FrameServer`.

    Example
    =======

    >>> frame_server = FrameServer(address="localhost", port=54321)
    >>> atoms = FaceCenteredCubic(directions=[[1, 0, 0], [0, 1, 0], [0, 0, 1]], symbol="Cu", size=(2, 2, 2), pbc=True)
    >>> dynamics = Langevin(atoms, timestep=0.5, temperature=300, friction=1.0)
    >>> dynamics.attach(send_ase_frame(atoms, frame_server), interval=2)
    """

    def send():
        frame = ase_to_frame_data(ase_atoms)
        frame_server.send_frame(send.frame_index, frame)
        send.frame_index = send.frame_index + 1

    send.frame_index = 0

    return send
