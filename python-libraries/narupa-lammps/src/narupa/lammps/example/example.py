import MDAnalysis
import numpy as np
from MDAnalysis.tests.datafiles import PSF, DCD  # test trajectory
from MDAnalysis import Universe
from MDAnalysis.topology.guessers import guess_atom_element

from narupa.protocol.trajectory import FrameData

element_index = {
    'H': 1,
    'C': 6,
    'N': 7,
    'O': 8,
    'S': 16
}

#U contains both topology and  positions
u = MDAnalysis.Universe(PSF, DCD)  # always start with a Universe

import time

from narupa.trajectory import FrameServer
from narupa.mdanalysis import mdanalysis_to_frame_data


frameServer = FrameServer(address='localhost', port=54321)
###Send topolgy once at the begining of the server
# Get topolgy in the right grpc format
#takes u because it likes things in mdalaysis format
topology_data = mdanalysis_to_frame_data(u, topology=True, positions=False)
#print(topology_data)

frame_index = 0
#Now actually send the frame, only contains the topology at this stage
frameServer.send_frame(frame_index, topology_data)

print("Starting Trajectory Server")
#frame_data = mdanalysis_to_frame_data(u, topology=False, positions=True)
#print(frame_data)


def lammps_to_frame_data(u: Universe, topology=True, positions=True) -> FrameData:
    frame_data = FrameData()

    if topology:
        for residue in u.residues:
            frame_data.arrays['residue.id'].string_values.values.append(residue.resname)
            frame_data.arrays['residue.chain'].index_values.values.append(residue.segment.ix)

        for atom in u.atoms:
            frame_data.arrays['atom.id'].string_values.values.append(atom.name)
            element = element_index[guess_atom_element(atom.name)]
            frame_data.arrays['atom.element'].index_values.values.append(element)
            frame_data.arrays['atom.residue'].index_values.values.append(atom.residue.ix)

        for bond in u.bonds:
            frame_data.arrays['bond'].index_values.values.append(bond.atoms[0].ix)
            frame_data.arrays['bond'].index_values.values.append(bond.atoms[1].ix)

    if positions:
        positions = np.multiply(0.1, np.ndarray.flatten(u.atoms.positions))
        frame_data.arrays["atom.position"].float_values.values.extend(positions)

    return frame_data



while True:

    for frame in u.trajectory:
        #Frame data is in grpc format
        frame_data = lammps_to_frame_data(u, topology=False, positions=True)
        print("FRAME STUFF",frame_index,frame_data)

        frameServer.send_frame(frame_index, frame_data)
        time.sleep(1.0 / 30.0)
        frame_index = frame_index + 1


class LammpsThing:
    def __init__(self):
        self.frame_server = FrameServer(address=..., port=...)
        self.frame_index = 0

    def hook(self, lmp):
        l = lammps(ptr=lmp)
        ...
        frame_data = lammps_to_frame_data(...)
        self.frame_server.send_frame(self.frame_index, frame_data)
        self.frame_index += 1
