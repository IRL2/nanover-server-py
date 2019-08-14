"""
LAMMPS python integration with Narupa
This program can be run as a standalone using dummy data or from within LAMMPS
using the python_invoke/fix command as demonstrated in the example LAMMPS inputs.
"""
import ctypes
from ctypes import c_int, c_double
import logging
from typing import Optional, Dict, Tuple, List

import numpy as np
try:
    from lammps import lammps
except ImportError:
    logging.info('lammps failed to import', exc_info=True)

from narupa.trajectory import FrameServer, FrameData
from narupa.trajectory.frame_data import POSITIONS, ELEMENTS

# IMD related imports
from narupa.imd.imd_force import calculate_imd_force
from narupa.imd.imd_server import ImdServer
from narupa.imd.imd_service import ImdService
from narupa.imd.particle_interaction import ParticleInteraction

# LAMMPS works with arbitrary masses, so we need to convert it to a nuclear number
# This list is a best guess for atom types, but won't work for isotopes for now.
ELEMENT_INDEX_MASS = {
    1:   1,
    3:   1,
    4:   2,
    7:   3,
    9:   4,
    11:  5,
    12:  6,
    14:  7,
    16:  8,
    19:  9,
    20:  10,
    23:  11,
    24:  12,
    27:  13,
    28:  14,
    31:  15,
    32:  16,
    35:  17,
    39:  19,
    40:  18,
    45:  21,
    48:  22,
    51:  23,
    52:  24,
    55:  25,
    56:  26,
    59:  27,
    64:  29,
    65:  30,
    70:  31,
    73:  32,
    75:  33,
    79:  34,
    80:  35,
    84:  36,
    85:  37,
    88:  38,
    89:  39,
}

# Check what units are being used in LAMMPS using this dict
# For now support converting lj and real, for full unit list see https://lammps.sandia.gov/doc/units.html
# List goes type, postion, force need to convert into nm and kj/mol/nm

LAMMPS_UNITS_CHECK = {
    # Lenard jones: Is unitless, everything is set to 1
    1.0           : ["lj", 1, 1],
    # Real:
    # Distance: 1 angstrom- > nm (10)
    # Force:    kj/mol/angstrom -> kcal/mol/nm (4.1840 *10) (Confirmed in MDanaysis)
    95.306976368  : ["real", 10, 41.840],
    # Metal:
    # Distance: angstrom -> nm, (10)
    # Force: eV/angstrom -> kcal/mol/nm (96.485*10) (Confirmed in MDanalysis)
    4.135667403e-3: ["metal", 10, 964.85],
    # SI:
    # Distance: meters ->nm (10^-9)
    # Force: Newtons q-> kcal/mol/nm (602214128999.9999)
    6.62606896e-34: ["si", 10**-9, 602214128999.9999],
    # cgs:
    # Distance: centemters -> nm
    # Froce: dyne (1/100000 Newtons) -> kj/mol*nm
    6.62606896e-27: ["cgs", 10**-7, 6022141.289999999],
    # Electron:
    # Distance: Bohr -> nm
    # Force: Hartree/Bohr (2625.50 / 0.0529117) ->  kj/mol*nm
    0.1519829846  : ["electron", 0.05529177, 49620.4053],
    # Mirco:
    # Distance: mircometers -> nm,
    # Force: pircogram-micrometer/microsecond^2 -> Newtons (1/1000000000000 *((1/1000000)/(1/1000000)^2) ->  kj/mol*nm
    6.62606896e-13: ["micro", 1000, 60221.41289999999],
    # Nano:
    # Distance: nanometers,
    # Force: atoogram-nanometer/nanosecond^2  -> Newtons (1/1e-12 *((1/1e-9)/(1/1e-9)^2) ->  kj/mol*nm
    6.62606896e-4:  ["nano", 1.0,602214128.9999999]
}



class DummyLammps:
    """
    A fake lammps object intended just for unit testing the lammps code
    without having to have lammps installed on a server
    """
    def __init__(self, n_atoms: int = None):
        _DEFAULT_ATOMS = 3
        self.n_atoms = n_atoms if n_atoms is not None else _DEFAULT_ATOMS

    def gather_atoms(self, array_type: str, dummy_variable, array_shape):
        """
        This routine generates fake ctypes to mimic lammps internal pointers

        :param array_type: determines the type of data that should be replicated
        :param dummy_variable: Unused here, only relevant to lammps
        :param array_shape: Unused here, only relevant to lammps
        :return: matrix data_array that contains all the dummy data
        """
        if array_type is "x":
            data_array = (ctypes.c_double * (3 * self.n_atoms))(*range(3 * self.n_atoms))
        elif array_type is "f":
            data_array = (ctypes.c_double * (3 * self.n_atoms))(*range(3 * self.n_atoms))
        elif array_type is "type":
            data_array = (ctypes.c_int * (self.n_atoms))(*range(1,1))
            # All atoms have the same type for DummyLammps testing
            for i in range(self.n_atoms):
                data_array[i] = 1
        else:
            logging.error('Unknown array type asked for in dummyLammps.gather_atoms')
            exit()

        return data_array

    def scatter_atoms(self, array_type, dummy_variable, array_shape, data_array):
        """
        This routine mimics L.scatter_atoms, in the dummy case it does nothing
        """

    def extract_global(self, types: str, number_type):
        '''
        Generate dummy element list for testing

        replicates L.extract_global("ntypes", 0)
        :param types: LAMMPS global variable that is needed
        :param number_type: unused
        :return:
        '''
        if types is "ntypes":
            dummy_element_list = 1
        elif types is "hplanck":
            dummy_element_list = 95.306976368
        else:
            logging.error('Unknown array type asked for in dummyLammps.extract_global')
            exit()
        return dummy_element_list

    def get_natoms(self):
        '''
        Generate dummy element list for testing
        '''
        return self.n_atoms

    def extract_atom(self, types: str, number_type):
        '''
        `Generate dummy element list for testing
        :param types: string that indicates the info that should be passed
        :param number_type: unused parameter to indicate integer float. etc
        :return: dummy_element_list
        '''
        if types is "mass":
            # initialise dummy type
            dummy_element_type = ctypes.c_double * 2
            # great array of dummy type
            dummy_element_list = dummy_element_type()
            # For some reason masses have a blank [0] value in LAMMPS
            dummy_element_list[0] = 0
            dummy_element_list[1] = 1
        else:
            logging.error('Unknown array type asked for in dummyLammps.extract_atom')
            exit()
        return dummy_element_list


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

    The translation of the LAMMPS c_type pointers into FrameData is done by  np.fromiter which
    allows a quick way of allocating a numpy array.

    The main lammps_hook routine will check if it is being run from within LAMMPS or as a
    stand alone program and determine if it should use dummy variables (manipulate_dummy_arrays)
    or ones available from within LAMMPS (manipulate_lammps_arrays).
    """

    def __init__(self, traj_port: int = 8080, imd_port: int = 8081, address: str = "[::]"):
        """
        Items that should be initialised on instantiation of lammpsHook class
        The MPI routines are essential to stop thread issues that cause internal
        LAMMPS crashes
        """
        # Start frame server, must come before MPI
        # TODO raise exception when this fails, i.e if port is blocked
        self.frame_server = FrameServer(address=address, port=traj_port)
        self.imd_server = ImdServer(address=address, port=imd_port)
        self.frame_index = 0
        self.frame_loop = 0

        try:
            self.frame_data = FrameData()
        except Exception as e:
            raise Exception("Failed to load FrameData", e)

        # Load MPI routines, has to be performed here.
        from mpi4py import MPI
        self.comm = MPI.COMM_WORLD
        me = self.comm.Get_rank()
        nprocs = self.comm.Get_size()

        logging.basicConfig(level=logging.INFO)
        logging.info("Lammpshook initialised for NarupaXR")
        logging.info("MPI rank %s", me)
        logging.info("MPI n processors %s ", nprocs)
        logging.info("Trajectory Port %s ", traj_port)
        logging.info("Interactive Port %s ", imd_port)
        # TODO make it so that the simulation waits on connect as an option

        # Set some variables that do not change during the simulation
        self.n_atoms = None
        self.units = None
        self.units_type = None
        self.force_factor = None
        self.distance_factor = None
        self.default_atoms = 10


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

    def manipulate_lammps_array(self, matrix_type: str, L):
        """
        Gather Matrix data from all LAMMPS MPI threads

        :param matrix_type: String identifying data to transmit, e.g x, v or f
        :param L: LAMMPS class that contains all the needed routines
        type :return: 3N matrix with all the data requested
        """
        data_array = L.gather_atoms(matrix_type, 1, 3)

        return data_array

    def gather_lammps_particle_types(self, L):
        '''
        Collect the particle list from LAMMPS, this may be atomistic or coarse grained
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
        atom_mass_type = [round(x) for x in atom_type_mass[0:ntypes+1]]
        # Convert to masses
        final_masses = [atom_mass_type[particle] for particle in atom_kind]
        final_masses = np.array(final_masses)
        # Convert to elements
        final_elements = [ELEMENT_INDEX_MASS.get(mass, 1) for mass in final_masses]
        return final_elements, final_masses

    def lammps_positions_to_frame_data(self,
                                       frame_data: FrameData,
                                       data_array: np.array) -> FrameData:
        """
        Convert the flat ctype.c_double data into the frame_data format. for the moment
        this assumes we are in LAMMPS real units. Its unclear at this stage if is possible
        to automatically detect this if the case is otherwise.
`
        :param data_array: Data to convert
        :param frame_data: frame data object
        """

        # Copy the ctype array to numpy for processing
        positions = np.fromiter(data_array, dtype=np.float, count=len(data_array))
        # Convert to nm
        positions = np.divide(positions, self.distance_factor)
        frame_data.arrays[POSITIONS] = positions

    @property
    def interactions(self) -> List[ParticleInteraction]:
        """
        Returns a shallow copy of the current interactions.
        This is copied from the ASE example, but reduces the abstracted degree
        """
        return self.imd_server.service.active_interactions

    def add_interaction_to_ctype(self, interaction_forces: np.array, lammps_forces):

        # initialise dummy type
        #ctype_3N_array = ctypes.c_double * (self.n_atoms * 3)
        # great array of dummy type
        #scatterable_array = ctype_3N_array()

        # Convert units fron narupa standard of  kj/mol/nm to whatever units LAMMPS is using
        # For real units types LAMMPS uses  Kcal/mole-Angstrom 4.14 for kj-> Kcal and 10x for nm -> Angstrom
        interaction_forces / self.force_factor
        # Flatten array into the ctype
        scatterable_array = interaction_forces.flatten()
        for idx in range(3*self.n_atoms):
            lammps_forces[idx] += scatterable_array[idx]

        return lammps_forces

    def return_array_to_lammps(self, matrix_type: str, scatter_array, L):
        L.scatter_atoms(matrix_type, 1, 3, scatter_array)

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
            try:
                L = DummyLammps(self.default_atoms)
            except Exception as e:
                # Many possible reasons for LAMMPS failures so for the moment catch all
                raise Exception("Failed to load DummyLammps", e)
        else:
            # Make sure LAMMPS object is callable
            try:
                L = lammps(ptr=lmp, comm=self.comm)
            except Exception as e:
                # Many possible reasons for LAMMPS failures so for the moment catch all
                raise Exception("Failed to load LAMMPS wrapper", e)

        # Check if we are in the first cycle of the MD to allocate static variables
        if self.n_atoms is None:
            self.n_atoms = L.get_natoms()
            # Use planks constant to work out what units we are working in
            self.units = L.extract_global("hplanck", 1)
            print("units", type(self.units),self.units)
            self.units_type = LAMMPS_UNITS_CHECK.get(self.units, None)[0]
            self.distance_factor = LAMMPS_UNITS_CHECK.get(self.units, None)[1]
            self.force_factor = LAMMPS_UNITS_CHECK.get(self.units, None)[2]
            print(self.units_type, self.force_factor, self.distance_factor)


        # Choose the matrix type that will be extracted
        positions = self.manipulate_lammps_array('x', L)
        # Copy the ctype array to numpy for processing
        positions3N = np.fromiter(positions, dtype=np.float, count=len(positions))
        # Convert to nm
        positions3N = np.multiply(self.distance_factor, positions3N).reshape(self.n_atoms, 3)

        # Collect forces from LAMMPS
        forces = self.manipulate_lammps_array('f', L)
        # Collect force vector from client
        interactions = self.interactions

        # get atom types and masses
        atom_type, masses = self.gather_lammps_particle_types(L)

        # Create numpy arrays with the forces to be added
        energy_kjmol, forces_kjmol = calculate_imd_force(positions3N, masses, interactions)
        # Add interactive force vector to lammps ctype
        forces = self.add_interaction_to_ctype(forces_kjmol, forces)

        # Return new force vector to LAMMPS
        self.return_array_to_lammps('f', forces, L)

        # Convert positions
        self.lammps_positions_to_frame_data(self.frame_data, positions)

        # Convert elements from list to frame data
        self.frame_data.arrays[ELEMENTS] = atom_type

        # Send frame data
        self.frame_server.send_frame(self.frame_index, self.frame_data)
        self.frame_index += 1

        # Print every 100 cycles if python interpreter is still running
        # This helps ensure that everything in lammps is continuing to run
        self.frame_loop += 1

        if self.frame_loop == 100:
            logging.info("LAMMPS python fix is running step %s", self.frame_index)
            #logging.info("FRAME STUFF %s %s", self.frame_index, self.frame_data.raw)
            self.frame_loop = 0
            print(self.frame_index)


