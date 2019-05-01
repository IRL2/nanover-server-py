"""
LAMMPS python integration with Naruapa
This program can be run as stand alone using dummy data or from within LAMMPS using
the python_invoke/fix command as demoed in the example LAMMPS inputs.
"""
import ctypes
import sys
import time
from ctypes import *
from pprint import pprint

import MDAnalysis
import mpi4py
import numpy as np
from lammps import lammps  # , PyLammps
from narupa.protocol.trajectory import FrameData
from narupa.trajectory import FrameServer
from numpy.core.multiarray import int_asbuffer

# Keep for converting internal LAMMPS atoms data into strings during testing
element_index = {
    'H': 1,
    'C': 6,
    'N': 7,
    'O': 8,
    'S': 16
}


class LammpsHook:
    def __init__(self):
        """
        Items that should be initialised on instantiation of lammpsHook class
        The MPI routines are essential to stop thread issues that cause internal
        LAMMPS crashes
        """

        # Load MPI routines
        from mpi4py import MPI
        self.comm = MPI.COMM_WORLD
        me = self.comm.Get_rank()
        nprocs = self.comm.Get_size()

        # Load frame server
        from narupa.protocol.trajectory import FrameData
        self.frame_server = FrameServer(address='localhost', port=54321)
        self.frame_index = 0
        print("Lammpshook initialised for NarupaXR")
        # TODO make it so that the simulation waits on connect as an option

    def testdebug(self):
        """
        Test routine to check correct python loading in LAMMPS
        keep for now, kill later

        :return: Nothing
        """

        try:
            L = lammps(ptr=lmp, comm=self.comm)
        except LammpsError as e:
            raise LammpsError("Failed to load LAMMPS wrapper", e)

        sys.exit(1)

        L = lammps(comm=comm)  # ptr=lmp, comm=comm)
        n_atoms = L.get_natoms()
        print("In class testy", "Atoms : ", n_atoms)

    def ManipulateLammpsArray(self, MatType, L):
        """
        Gather Matrix data from all LAMMPS MPI threads

        :param MatType: String identifying data to transmit, e.g x, v or f
        :param L: LAMMPS class that contains all the needed routines
        :return: 3N matrix called v with all the data requested
        """

        # n_local = L.extract_global('nlocal', 0)  # L.get_nlocal()
        # Hard to tell if LAMMPS python interpreter is working so for now print every step
        print("In LAMMPS array")
        n_atoms = L.get_natoms()
        v = L.gather_atoms(MatType, 1, 3)

        # This test case slowly translates the molecular system
        for idx in range(n_atoms):
            v[3 * idx + 0] += 0.0001000
            v[3 * idx + 1] *= 1.0000000
            v[3 * idx + 2] *= 1.0000000
        L.scatter_atoms(MatType, 1, 3, v)
        return v

    def ManipulateDummyArray(self, MatType):
        """
        This routine mimics LAMMPS cytpes for easy debugging
        Generate dummy ctype double array of 3N particles
        TODO convert this to a full dummy LAMMPS class

        :param MatType: For the moment doesnt do anything
        :return: 3N matrix v that contains all the dummy data
        """
        n_atoms = 10
        v = (ctypes.c_double * (3 * n_atoms))(*range(3 * n_atoms))
        print(v[1], v[2], v[3])
        return v

    # def LammpsFrameDataArray(self):
    #     # Convert lammps array to gprc frame format.
    #     print()
    #
    # def LammpsFrameDataValue(self):
    #     # Convert lammps value to gprc frame format.
    #     print()

    def lammps_to_frame_data(self, v, topology=True, positions=True) -> FrameData:
        """
        Convert the flat ctype.c_double data into the framedata format.

        :param v: Data to convert
        :param topology: Check if data is topolgical
        :param positions: Check if data is positional
        :return: overwrite data in v matrix with new formatted framedata
        """
        try:
            frame_data = FrameData()
        except Exception as e:
            raise Exception("Failed to load framedata", e)
            sys.exit(1)
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
            # print("Ndarray",np.ndarray(v))

            # positions = Array.from_address(ctypes.addressof(v.contents))
            # positions = np.multiply(0.1, np.frombuffer(v))
            # Copy the ctype array to numpy for processing
            positions = np.fromiter(v, dtype=np.float, count=len(v))
            # Convert to nm
            positions = np.multiply(0.1, positions)
            frame_data.arrays["atom.position"].float_values.values.extend(positions)

        return frame_data

    # LammpsHook passes data between python and the Lammps binary
    def LammpsHook(self, lmp=None):
        """
        LammpsHook is the main routine that is run within LAMMPS MD
        steps. It checks that LAMMPS python wrapper is callable
        and then attempts to extract a 3N matrix of atomic data

        :param lmp: LAMMPS object data, only populated when running from within LAMMPS
        """

        # Checks if LAMMPS variable is being passed
        # If not assume we are in interactive mode
        if lmp is None:
            print("Running without lammps, assuming interactive debugging")
        else:
            # Make sure LAMMPS object is callable
            try:
                L = lammps(ptr=lmp, comm=self.comm)
            except LammpsError as e:
                # Many reasons for LAMMPS failures so for the moment catch all
                raise LammpsError("Failed to load LAMMPS wrapper", e)
                sys.exit(1)

        # mass = L.extract_atom("mass",2)
        # xp = L.extract_atom("x",3)
        # print("Natoms, mass, x[0][0] coord =",n_atoms,mass[1],xp[0][0])
        # temp = L.extract_compute("thermo_temp",0,0)
        # print("Temperature from compute =",temp)

        # Choose the matrix type that will be extracted
        MatType = "x"
        # If not in LAMMPS run dummy routine
        if lmp is None:
            v = self.ManipulateDummyArray(MatType)
        else:
            v = self.ManipulateLammpsArray(MatType, L)

        self.frame_data = self.lammps_to_frame_data(v, positions=True, topology=False)

        # print("FRAME STUFF \n", self.frame_index, "\n", self.frame_data)
        self.frame_server.send_frame(self.frame_index, self.frame_data)
        self.frame_index += 1

        # Scatter data back to lammps processors
        # if lmp is not None:
        # L.scatter_atoms(MatType,1,3,v)


# Test call of the routine when running outside of lammps
def main():
    H = LammpsHook()
    print("Starting Trajectory Server")
    # frameServer = FrameServer(address='localhost', port=54321)
    while True:
        H.LammpsHook()
        print("FRAME STUFF", H.frame_index, H.frame_data)
        time.sleep(1.0 / 10.0)


if __name__ == '__main__':
    main()
