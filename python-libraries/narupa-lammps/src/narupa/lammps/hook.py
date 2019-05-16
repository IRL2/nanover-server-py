"""
LAMMPS python integration with Narupa
This program can be run as a standalone using dummy data or from within LAMMPS
using the python_invoke/fix command as demonstrated in the example LAMMPS inputs.
"""
import ctypes
import logging

import numpy as np
from lammps import lammps  # , PyLammps

from narupa.protocol.trajectory import FrameData
from narupa.trajectory import FrameServer, FrameData
from narupa.trajectory.frame_data import POSITIONS, ELEMENTS

# Keep for converting internal LAMMPS atoms data into strings during testing
element_index_mass = {
    1 :'H',
    12:'C',
    14:'N',
    16:'O',
    32:'S'
}


def manipulate_dummy_array(matrix_type, n_atoms):
    """
    This routine mimics LAMMPS cytpes for easy debugging
    Generate dummy ctype double array of 3N particles
    TODO convert this to a full dummy LAMMPS class

    :param matrix_type: For the moment doesnt do anything
    :param n_atoms: Number of atoms detemrines dimension of array
    :return: 3N matrix data_array that contains all the dummy data
    """
    data_array = (ctypes.c_double * (3 * n_atoms))(*range(3 * n_atoms))
    print(data_array[1], data_array[2], data_array[3])
    return data_array


class LammpsHook:
    """
    lammps_hook is a series of routines the can communicate with the LAMMPS program through
    its python interpreter. Upon initialisation, MPI is set up along with the frame server.
    The LAMMPS data is collected across all processors using GATHER and SCATTER routines
    that require mpi4py to respect the internal processor rank of LAMMPS.

    The variables that can currently be accessed are
    x : positions
    v : velocities
    f : forces

    We hook into the LAMMPS MD loop just after the forces are calculated, however setting all
    forces to zero causes the energy in LAMMPS to become large. To check the effect of this code
    on the internal state of LAMMPS we recommend trying to set the velocities to zero, effectively
    freezing the system. In the case below we are slowly translating all the atoms in the system
    along the x direction.

    The translation of the LAMMPS c_type pointers into framedata is done by  np.fromiter which
    allows a quick way of allocating a numpy array.

    The main lammps_hook routine will check if it is being run from within LAMMPS or as a
    stand alone program and determine if it should use dummy variables (manipulate_dummy_arrays)
    or ones available from within LAMMPS (manipulate_lammps_arrays).
    """

    def __init__(self):
        """
        Items that should be initialised on instantiation of lammpsHook class
        The MPI routines are essential to stop thread issues that cause internal
        LAMMPS crashes
        """
        # Start frame server, must come before MPI
        port_no = 54321
        self.frame_server = FrameServer(address='[::]', port=port_no)
        self.frame_index = 0
        self.frame_loop = 0

        # Load MPI routines, has to be performed here.
        from mpi4py import MPI
        self.comm = MPI.COMM_WORLD
        me = self.comm.Get_rank()
        nprocs = self.comm.Get_size()

        logging.basicConfig(level=logging.INFO)
        logging.info("Lammpshook initialised for NarupaXR")
        logging.info("MPI rank %s", me)
        logging.info("MPI n processors %s ", nprocs)
        logging.info("Port %s ", port_no)
        # TODO make it so that the simulation waits on connect as an option

    def test_debug(self):
        """
        Test routine to check correct python loading in LAMMPS
        TODO remove once more robust testing of loading in LAMMPS is developed
        """
        try:
            L = lammps(ptr=lmp, comm=self.comm)
        except Exception as e:
            raise Exception("Failed to load LAMMPS wrapper", e)

        n_atoms = L.get_natoms()
        print("In class testy", "Atoms : ", n_atoms)

    def manipulate_lammps_array(self, matrix_type, L):
        """
        Gather Matrix data from all LAMMPS MPI threads

        :param matrix_type: String identifying data to transmit, e.g x, v or f
        :param L: LAMMPS class that contains all the needed routines
 type :return: 3N matrix with all the data requested
        """

        # n_local = L.extract_global('nlocal', 0)  # L.get_nlocal()
        # Hard to tell if LAMMPS python interpreter is working so for now print every step
        n_atoms = L.get_natoms()
        data_array = L.gather_atoms(matrix_type, 1, 3)

        # This test case slowly translates the molecular system
        for idx in range(n_atoms):
            data_array[3 * idx + 0] += 0.0001000
            data_array[3 * idx + 1] *= 1.0000000
            data_array[3 * idx + 2] *= 1.0000000
        L.scatter_atoms(matrix_type, 1, 3, data_array)
        return data_array

    def gather_lammps_particle_types(self, L):
        '''
        Collect the paticle list from LAMMPS, this may be atomistic or coarse grained
        particles by looking up the ID of the atoms and that ID's corresponding mass
        from the atoms_elements dict.

        :param L: LAMMPS class that contains all the needed routines
        :return: 1N matrix with all the data requested
        '''

        # Extract the number of atoms types in the system.
        ntypes = L.extract_global("ntypes", 0)
        # Extract the masses of the types, 1D float of home many
        # mass types were defined in input. Indexed from 1 not zero in lammps
        atom_type_mass = L.extract_atom("mass", 2)
        # Gather the atom types, 1D int of n atoms length.
        atom_kind = L.gather_atoms("type", 0, 1)


        # Atom mass is indexed from 1 in lammps for some reason.
        # Create a new list rounded to the nearest mass integer
        atom_mass_type = list(atom_type_mass[1:ntypes+1])
        atom_mass_type = [round(x) for x in atom_mass_type]

        atom_elements =[atom_mass_type[particle-1] for particle in atom_kind]


        final_elements =[element_index_mass[mass] for mass in atom_elements]
        #print(atom_mass_type)
        #print(atom_kind[0:10])
        #print(final_elements[0:10])
        return final_elements

    def lammps_array_to_frame_data(self, data_array, frame_data) -> FrameData:
        """
        Convert the flat ctype.c_double data into the framedata format.
`
        :param data_array: Data to convert
        :param topology: Check if data is topological
        :param positions: Check if data is positional
        :return: overwrite data in data_array matrix with new formatted framedata
        """

        # Copy the ctype array to numpy for processing
        positions = np.fromiter(data_array, dtype=np.float, count=len(data_array))
        # Convert to nm
        positions = np.multiply(0.1, positions)
        frame_data.arrays[POSITIONS] = positions

#    def lammps_elements_to_frame_data(self, data_array,  frame_data) -> FrameData:
#        """
#        Convert the flat ctype.c_double data into the framedata format.
#
#        :param data_array: Data to convert
#        :param topology: Check if data is topological
#        :param positions: Check if data is positional
#        :return: overwrite data in data_array matrix with new formatted framedata
#        """
#
#
        return frame_data


    # lammps_hook passes data between python and the Lammps binary
    def lammps_hook(self, lmp=None):
        """
        lammps_hook is the main routine that is run within LAMMPS MD
        steps. It checks that the LAMMPS python wrapper is callable
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
            except Exception as e:
                # Many possible reasons for LAMMPS failures so for the moment catch all
                raise Exception("Failed to load LAMMPS wrapper", e)

        # Choose the matrix type that will be extracted
        matrix_type = "x"
        n_atoms_dummy = 10
        # If not in LAMMPS run dummy routine
        if lmp is None:
            data_array = manipulate_dummy_array(matrix_type, n_atoms_dummy)
        else:
            data_array = self.manipulate_lammps_array(matrix_type, L)
            atom_type = self.gather_lammps_particle_types(L)

        # Convert collected data to frame data
        try:
            frame_data = FrameData()
        except Exception as e:
            raise Exception("Failed to load framedata", e)

        # Convert positions
        self.frame_data = self.lammps_array_to_frame_data(data_array, frame_data)

        # Remove once new pretend lammps object is made
        # Convert elements
        if lmp is not None:
            self.frame_data.arrays[ELEMENTS] = atom_type
            #self.frame_data = self.lammps_elements_to_frame_data(atom_type, frame_data)

        # Send frame data
        # print("FRAME STUFF \n", self.frame_index, "\n", self.frame_data)
        self.frame_server.send_frame(self.frame_index, self.frame_data)
        self.frame_index += 1

        self.frame_loop += 1
        # Print every 100 cycles if python interperater is still running
        if self.frame_loop == 10 :
            logging.info("LAMMPS python fix is running step %s", self.frame_index)
            logging.info("FRAME STUFF %s %s", self.frame_index, self.frame_data.raw)
            self.frame_loop = 0

