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
import ctypes
from numpy.core.multiarray import int_asbuffer

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


class LammpsHook:
    def __init__(self):
        #Load MPI routines
        from mpi4py import MPI
        self.comm = MPI.COMM_WORLD
        me = self.comm.Get_rank()
        nprocs = self.comm.Get_size()
        #Load frame server
        from narupa.protocol.trajectory import FrameData
        self.frame_server = FrameServer(address='localhost', port=54321)
        self.frame_index = 0
        print("Lammpshook initalised")

    # Test routine to check correct loading, keep for now, kill later
    def testy(self):
        try:
            L = lammps(ptr=lmp, comm=self.comm)
        except LammmpsError:
            print("Failed to load lammps wrapper")
        L = lammps(comm=comm)# ptr=lmp, comm=comm)
        print("inclass.testy")
        n_atoms = L.get_natoms()
        print("Atoms : ", n_atoms)

    def ManipulateLammpsArray(self, MatType, L):
        #n_local = L.extract_global('nlocal', 0)  # L.get_nlocal()
        print("in LAMMPS array")
        n_atoms = L.get_natoms()
        v = L.gather_atoms(MatType, 1, 3)
        for idx in range(n_atoms):
            v[3 * idx + 0] *= 0.0000001
            v[3 * idx + 1] *= 0.0000001
            v[3 * idx + 2] *= 0.0000001
        return v

    # This routine mimics LAMMPS cytpes for easy debugging
    def ManipulateDummyArray(self, MatType):
        import ctypes
        n_atoms = 10
        v = (ctypes.c_double*(3*n_atoms))(*range(3*n_atoms))
        print(v[1],v[2],v[3])
        return v

    # def LammpsFrameDataArray(self):
    #     # Convert lammps array to gprc frame format.
    #     print()
    #
    # def LammpsFrameDataValue(self):
    #     # Convert lammps value to gprc frame format.
    #     print()

    def lammps_to_frame_data(self, v, topology=True, positions=True) -> FrameData:
        try:
            frame_data = FrameData()
        except Exception as e:
            print("Failed to load framedata")
            print(e)


        # if topology:
        #     for residue in u.residues:
        #         frame_data.arrays['residue.id'].string_values.values.append(residue.resname)
        #         frame_data.arrays['residue.chain'].index_values.values.append(residue.segment.ix)
        #
        #     for atom in u.atoms:
        #         frame_data.arrays['atom.id'].string_values.values.append(atom.name)
        #         element = element_index[guess_atom_element(atom.name)]
        #         frame_data.arrays['atom.element'].index_values.values.append(element)
        #         frame_data.arrays['atom.residue'].index_values.values.append(atom.residue.ix)
        #
        #     for bond in u.bonds:
        #         frame_data.arrays['bond'].index_values.values.append(bond.atoms[0].ix)
        #         frame_data.arrays['bond'].index_values.values.append(bond.atoms[1].ix)

        if positions:
            #print("Ndarray",np.ndarray(v))

            #positions = Array.from_address(ctypes.addressof(v.contents))
            #positions = np.multiply(0.1, np.frombuffer(v))
            # Copy the ctype array to numpy for processing
            positions = np.fromiter(v, dtype=np.float, count=len(v))
            # Convert to nm
            positions = np.multiply(0.1, positions)
            frame_data.arrays["atom.position"].float_values.values.extend(positions)

        return frame_data

    # LammpsHook passes data between python and the Lammps binary
    def LammpsHook(self, lmp=None):
        # Check if LAMMPS variable is being passed
        # If not assume we are in interactive mode
        if lmp is None:
            print("Running without lammps, assuming interactive")
        else:
            # Make sure lammps object is callable
            try:
                L = lammps(ptr=lmp, comm=self.comm)
            except Exception as e:
                print("Failed to load LAMMPS wrapper")
                print(e)

        # mass = L.extract_atom("mass",2)
        # xp = L.extract_atom("x",3)
        # print("Natoms, mass, x[0][0] coord =",n_atoms,mass[1],xp[0][0])
        # temp = L.extract_compute("thermo_temp",0,0)
        # print("Temperature from compute =",temp)

        # n3=3*n_atoms
        MatType= "v"
        if lmp is None:
            v = self.ManipulateDummyArray(MatType)
        else:
            v = self.ManipulateLammpsArray(MatType, L)
        self.frame_data = self.lammps_to_frame_data(v, positions=True, topology=False)
        #print("FRAME STUFF \n", self.frame_index, "\n", self.frame_data)
        if lmp is not None:
            self.frame_server.send_frame(self.frame_index, self.frame_data)
        self.frame_index += 1


#Test call of the routine
#frameServer = FrameServer(address='localhost', port=54321)
#H = LammpsHook()
#while True:
#for x in range(0,10):
#    H.LammpsHook()
#
#    # Frame data is in grpc format
#    print("FRAME STUFF", H.frame_index, H.frame_data)
#    frameServer.send_frame(H.frame_index, H.frame_data)
#    time.sleep(1.0 / 30.0)
#    #frame_index = frame_index + 1
