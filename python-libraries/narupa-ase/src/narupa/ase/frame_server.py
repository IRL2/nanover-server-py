# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.

"""
Methods for transmitting a simulation frame from ASE.
"""

from ase import Atoms
from narupa.trajectory import FrameServer

from . import ase_to_frame_data


def send_ase_frame(ase_atoms: Atoms, frame_server: FrameServer):
    """
    Hook to transmit the current state of an ASE Atoms as a frame.
    :param ase_atoms: ASE atoms
    :param frame_server: The frame server on which to send the produced frame.
    """
    def send():
        frame = ase_to_frame_data(ase_atoms)
        frame_server.send_frame(send.frame_index, frame)
        send.frame_index = send.frame_index + 1

    send.frame_index = 0

    return send
