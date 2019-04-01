from narupy.ase.frame import ase_atoms_to_frame_data
from narupy.ase.topology import ase_atoms_to_topology_data
from ase import Atoms

def NarupaASE(ase_atoms : Atoms, frameServer):
    def send():
        frame = ase_atoms_to_frame_data(ase_atoms)
        topology = ase_atoms_to_topology_data(ase_atoms)
        frameServer.send_topology(send.frame_index, topology)
        frameServer.send_frame(send.frame_index, frame)
        send.frame_index = send.frame_index + 1
    send.frame_index = 0

    return send

