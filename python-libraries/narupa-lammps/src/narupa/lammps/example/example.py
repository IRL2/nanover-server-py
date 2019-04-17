import MDAnalysis
import numpy as np
from MDAnalysis.tests.datafiles import PSF, DCD  # test trajectory
from MDAnalysis import Universe
from MDAnalysis.topology.guessers import guess_atom_element
# from __future__ import print_function
from lammps import lammps  # , PyLammps
import mpi4py
# from lammps import lammps #, PyLammps
import sys
from pprint import pprint
from ctypes import *

from narupa.protocol.trajectory import FrameData

element_index = {
    'H': 1,
    'C': 6,
    'N': 7,
    'O': 8,
    'S': 16
}

# U contains both topology and  positions
# u = MDAnalysis.Universe(PSF, DCD)  # always start with a Universe
#
import time
#
from narupa.trajectory import FrameServer
# from narupa.mdanalysis import mdanalysis_to_frame_data
#
# frameServer = FrameServer(address='localhost', port=54321)
# ###Send topolgy once at the begining of the server
# # Get topolgy in the right grpc format
# # takes u because it likes things in mdalaysis format
# topology_data = mdanalysis_to_frame_data(u, topology=True, positions=False)
# # print(topology_data)
#
# frame_index = 0
# # Now actually send the frame, only contains the topology at this stage
# frameServer.send_frame(frame_index, topology_data)

print("Starting Trajectory Server")


# frame_data = mdanalysis_to_frame_data(u, topology=False, positions=True)
# print(frame_data)


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

#L = lammps(ptr=lmp, comm=comm)

class LammpsHook:
    def __init__(self):
        #Load MPI routines
        from mpi4py import MPI
        self.comm = MPI.COMM_WORLD
        me = self.comm.Get_rank()
        nprocs = self.comm.Get_size()
        #Load frame server
        #self.frame_server = FrameServer(address='localhost', port=54321)
        #self.frame_index = 0

    def testy(self):
        #from mpi4py import MPI
        #comm = MPI.COMM_WORLD
        #me = comm.Get_rank()
        #nprocs = comm.Get_size()
        try:
            L = lammps(ptr=lmp, comm=self.comm)
        except LammmpsError:
            print("Failed to load lammps wrapper")
        L = lammps(comm=comm)# ptr=lmp, comm=comm)
        print("inclass.testy")
        n_atoms = L.get_natoms()
        print("Atoms : ", n_atoms)

    def LammpsHook(self, lmp):
        # Make sure lammps object is callable
        try:
            L = lammps(ptr=lmp, comm=self.comm)
        except Exception as e:
            print("Failed to load lammps wrapper")
            print(e)
        n_atoms = L.get_natoms()
        # mass = L.extract_atom("mass",2)
        # xp = L.extract_atom("x",3)
        # print("Natoms, mass, x[0][0] coord =",n_atoms,mass[1],xp[0][0])
        # temp = L.extract_compute("thermo_temp",0,0)
        # print("Temperature from compute =",temp)

        n_local = L.extract_global('nlocal', 0)  # L.get_nlocal()
        # print(tmp_force.__dict__)
        # n3=3*n_atoms
        v = L.gather_atoms("v", 1, 3)
        for idx in range(n_atoms):
            # print(type(v[idx]))
            v[3 * idx + 0] *= 0.0000001
            v[3 * idx + 1] *= 0.0000001
            v[3 * idx + 2] *= 0.0000001
        L.scatter_atoms("v", 1, 3, v)

        ...
        # frame_data = lammps_to_frame_data(...)
        # self.frame_server.send_frame(self.frame_index, frame_data)
        # self.frame_index += 1


# while True:
#
#     for frame in u.trajectory:
#         # Frame data is in grpc format
#         frame_data = lammps_to_frame_data(u, topology=False, positions=True)
#         print("FRAME STUFF", frame_index, frame_data)
#
#         frameServer.send_frame(frame_index, frame_data)
#         time.sleep(1.0 / 30.0)
#         frame_index = frame_index + 1
