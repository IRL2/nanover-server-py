import MDAnalysis
from MDAnalysis.tests.datafiles import PSF, DCD  # test trajectory

u = MDAnalysis.Universe(PSF, DCD)  # always start with a Universe

from narupa.protocol.trajectory.frame_pb2 import FrameData
from narupa.protocol.topology.topology_pb2 import TopologyData

import time

import numpy as np

from narupy.trajectory.frame_server import FrameServer

element_index = {
    'H': 1,
    'C': 6,
    'N': 7,
    'O': 8,
    'S': 16
}

frameServer = FrameServer(address='localhost', port=54321)

topology_data = TopologyData()

for residue in u.residues:
    topology_data.arrays['residue.id'].string_values.values.append(residue.resname)
    topology_data.arrays['residue.chain'].index_values.values.append(residue.segment.ix)

for atom in u.atoms:
    topology_data.arrays['atom.id'].string_values.values.append(atom.name)
    element = element_index[MDAnalysis.topology.guessers.guess_atom_element(atom.name)]
    topology_data.arrays['atom.element'].index_values.values.append(element)
    topology_data.arrays['atom.residue'].index_values.values.append(atom.residue.ix)

for bond in u.bonds:
    topology_data.arrays['bond'].index_values.values.append(bond.atoms[0].ix)
    topology_data.arrays['bond'].index_values.values.append(bond.atoms[1].ix)

frame_index = 0
frameServer.send_topology(frame_index, topology_data)

print("Starting Trajectory Server")

while True:

    for frame in u.trajectory:
        frame_data = FrameData()

        positions = np.multiply(0.1, np.ndarray.flatten(u.atoms.positions))
        frame_data.arrays["atom.position"].float_values.values.extend(positions)

        frameServer.send_frame(frame_index, frame_data)
        time.sleep(1.0 / 30.0)
        frame_index = frame_index + 1
