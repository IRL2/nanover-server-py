# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.
"""
LAMMPS python integration with Narupa
This program can be run as a standalone using dummy data or from within LAMMPS
using the python_invoke/fix command as demonstrated in the example LAMMPS inputs.
"""
import ctypes
import functools
import logging
from typing import List
import numpy as np

try:
    from lammps import lammps
except ImportError:
    logging.info('lammps failed to import', exc_info=True)

from narupa.trajectory import FrameServer, FrameData
from narupa.trajectory.frame_data import PARTICLE_POSITIONS, PARTICLE_ELEMENTS

# IMD related imports
from narupa.imd.imd_force import calculate_imd_force
from narupa.imd.imd_server import ImdServer
from narupa.imd.particle_interaction import ParticleInteraction

# LAMMPS works with arbitrary masses, so we need to convert it to a nuclear number
# This list is a best guess for atom types, but won't work for isotopes for now.
# see issue 82 https://gitlab.com/intangiblerealities/narupa-protocol/issues/82
ELEMENT_INDEX_MASS = {
    1: 1,
    3: 1,
    4: 2,
    7: 3,
    9: 4,
    11: 5,
    12: 6,
    14: 7,
    16: 8,
    19: 9,
    20: 10,
    23: 11,
    24: 12,
    27: 13,
    28: 14,
    31: 15,
    32: 16,
    35: 17,
    39: 19,
    40: 18,
    45: 21,
    48: 22,
    51: 23,
    52: 24,
    55: 25,
    56: 26,
    59: 27,
    64: 29,
    65: 30,
    70: 31,
    73: 32,
    75: 33,
    79: 34,
    80: 35,
    84: 36,
    85: 37,
    88: 38,
    89: 39,
}

# Check what units are being used in LAMMPS using this dict
# For now support converting lj and real, for full unit list
# see https://lammps.sandia.gov/doc/units.html
# List goes type, postion, force need to convert into nm and kj/mol/nm

LAMMPS_UNITS_CHECK = {
    # Lenard jones: Is unitless, everything is set to 1
    0: ["lj", 1, 1],
    # Real:
    # Distance: 1 angstrom- > nm (10)
    # Force:    kj/mol/angstrom -> kcal/mol/nm (4.1840 *10) (Confirmed in MDanaysis)
    1: ["real", 10, 41.840],
    # Metal:
    # Distance: angstrom -> nm, (10)
    # Force: eV/angstrom -> kcal/mol/nm (96.485*10) (Confirmed in MDanalysis)
    2: ["metal", 10, 964.85],
    # SI:
    # Distance: meters ->nm (10^-9)
    # Force: Newtons q-> kcal/mol/nm (602214128999.9999)
    3: ["si", 10 ** -9, 602214128999.9999],
    # cgs:
    # Distance: centemters -> nm
    # Froce: dyne (1/100000 Newtons) -> kj/mol*nm
    4: ["cgs", 10 ** -7, 6022141.289999999],
    # Electron:
    # Distance: Bohr -> nm
    # Force: Hartree/Bohr (2625.50 / 0.0529117) ->  kj/mol*nm
    5: ["electron", 0.05529177, 49620.4053],
    # Mirco:
    # Distance: mircometers -> nm,
    # Force: pircogram-micrometer/microsecond^2 -> Newtons
    # (1/1000000000000 *((1/1000000)/(1/1000000)^2) ->  kj/mol*nm
    6: ["micro", 1000, 60221.41289999999],
    # Nano:
    # Distance: nanometers,
    # Force: atoogram-nanometer/nanosecond^2  -> Newtons
    # (1/1e-12 *((1/1e-9)/(1/1e-9)^2) ->  kj/mol*nm
    7: ["nano", 1.0, 602214128.9999999]
}
# store plank values as a list so that we don't do float lookups in a dict.
PLANK_VALUES = (
    1.0,
    95.306976368,
    4.135667403e-3,
    6.62606896e-34,
    6.62606896e-27,
    0.1519829846,
    6.62606896e-13,
    6.62606896e-4
)


class DummyLammps:
    """
    A fake lammps object intended just for unit testing the lammps code
    without having to have lammps installed on a server
    """

    def __init__(self, n_atoms_in_dummy: int = None):
        # Set a default atom length for tests
        _DEFAULT_ATOMS = 3
        self.n_atoms = n_atoms_in_dummy if n_atoms_in_dummy is not None else _DEFAULT_ATOMS

    def gather_atoms(self, array_type: str, _dummy_variable, _array_shape):
        """
        This routine generates fake ctypes to mimic lammps internal pointers

        :param array_type: determines the type of data that should be replicated
        :param _dummy_variable: Unused here, only relevant to lammps
        :param _array_shape: Unused here, only relevant to lammps
        :return: matrix data_array that contains all the dummy data
        """
        empty_list = []
        if array_type == "x":
            data_array = (ctypes.c_double * (3 * self.n_atoms))(*range(3 * self.n_atoms))
        elif array_type == "f":
            data_array = (ctypes.c_double * (3 * self.n_atoms))(*empty_list)
        elif array_type == "type":
            data_array = (ctypes.c_int * self.n_atoms)(*empty_list)
            # All atoms have the same type for DummyLammps testing
            for i in range(self.n_atoms):
                data_array[i] = 1
        else:
            raise Exception('Unknown array type asked for in dummyLammps.gather_atoms')

        return data_array

    def scatter_atoms(self, _array_type, _dummy_variable, _array_shape, __data_array):
        """
        This routine mimics lammp_class.scatter_atoms, in the dummy case it does nothing
        """
        return None

    def extract_global(self, types: str, _number_type):
        """
        Generate dummy element list for testing

        replicates lammp_class.extract_global("ntypes", 0)
        :param types: LAMMPS global variable that is needed
        :param _number_type: unused
        :return:
        """
        if types == "ntypes":
            dummy_element_list = 1
        elif types == "hplanck":
            dummy_element_list = 95.306976368
        else:
            raise Exception('Unknown array type asked for in dummyLammps.extract_global')

        return dummy_element_list

    def get_natoms(self):
        """
        Generate dummy element list for testing
        """
        return self.n_atoms

    def extract_atom(self, types: str, _number_type):
        """
        Generate dummy element list for testing
        :param types: string that indicates the info that should be passed
        :param _number_type: unused parameter to indicate integer float. etc
        :return: dummy_element_list
        """
        if types == "mass":
            # initialise dummy type
            dummy_element_type = ctypes.c_double * 2
            # great array of dummy type
            dummy_element_list = dummy_element_type()
            # For some reason masses have a blank [0] value in LAMMPS
            dummy_element_list[0] = 0
            dummy_element_list[1] = 1
        else:
            raise Exception('Unknown array type asked for in dummyLammps.extract_atom')

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
    topology_loop = True

    def __init__(self, traj_port: int = None, imd_port: int = None, address: str = "[::]"):
        """
        Items that should be initialised on instantiation of lammpsHook class
        The MPI routines are essential to stop thread issues that cause internal
        LAMMPS crashes
        """
        logging.basicConfig(level=logging.INFO)
        me = 0
        try:
            from mpi4py import MPI
            self.comm = MPI.COMM_WORLD
            me = MPI.COMM_WORLD.Get_rank()
            nprocs = MPI.COMM_WORLD.Get_size()
            self.me = me
            self.nprocs = nprocs

            if me == 0:
                logging.info("MPI rank %s", me)
                logging.info("MPI n processors %s ", nprocs)
        except ImportError as err:
            logging.info("Didn't find mpi4py %s", err)
        # Start frame server, must come before MPI
        if me == 0:
            # TODO raise exception when this fails, i.e if port is blocked
            self.frame_server = FrameServer(address=address, port=traj_port)
            self.imd_server = ImdServer(address=address, port=imd_port)
            self.frame_index = 0
            self.frame_loop = 0

            try:
                self.frame_data = FrameData()
            except Exception as err:
                raise Exception("Failed to load FrameData", err)

            logging.info("Lammpshook initialised for NarupaXR")
            # Load MPI routines, has to be performed here.

            logging.info("Trajectory Port %s ", traj_port)
            logging.info("Interactive Port %s ", imd_port)

            # Set some variables that do not change during the simulation
            self.n_atoms = None
            self.units = None
            self.units_type = None
            self.force_factor = None
            self.distance_factor = None
            self.masses = None
            self.atom_type = None
            self.n_atoms_in_dummy = 10
            self.loop = 0

    def _try_or_except(func):
        """
        Function creates an except or try for various functions to overcome the issue
        of the LAMMPS interpreter not giving error messages correctly when an error is
        enocuntered.
        :param func: function to be decorated with a try or except statement
        :return: the original function but decorated
        """

        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            try:
                func(*args, **kwargs)
            except Exception as e:
                # Note args[0] is used to get around the issue of passing self to a decorator
                logging.info("Exception raised in calling function on proc %s ", args[0].me)
                logging.info("Exception thrown %s ", e)
            return func(*args, **kwargs)

        return wrapper

    def close(self):
        """
        Close ports to prevent blocking
        """
        logging.info("Closing Narupa server")
        self.frame_server.close()
        self.imd_server.close()

    @_try_or_except
    def manipulate_lammps_array(self, matrix_type: str, lammps_class):
        """
        Gather Matrix data from all LAMMPS MPI threads

        :param matrix_type: String identifying data to transmit, e.g x, v or f
        :param lammps_class: LAMMPS class that contains all the needed routines
        type :return: 3N matrix with all the data requested
        """

        data_array = lammps_class.gather_atoms(matrix_type, 1, 3)

        return data_array

    def gather_lammps_particle_types(self, lammps_class):
        """
        Collect the particle list from LAMMPS, this may be atomistic or coarse grained
        particles by looking up the ID of the atoms and that ID's corresponding mass
        from the atoms_elements dict.

        :param lammps_class: LAMMPS class that contains all the needed routines
        :return: 1N matrix with all the data requested
        """
        print("AAAAHAHHAHAHAHHAHAHAHAHAHAHHAHAHAH  ", self.me)
        # me = self.comm.Get_rank()
        # nprocs = self.comm.Get_size()
        # logging.info("gather MPI rank %s", me)
        # logging.info("gather MPI n processors %s ", nprocs)
        # Extract the number of atoms types in the system.
        ntypes = lammps_class.extract_global("ntypes", 0)
        logging.info("test 1 %s %s", ntypes, self.me)

        # Extract the masses of the types, 1D float of home many
        # mass types were defined in input. Indexed from 1 not zero in lammps
        atom_type_mass = lammps_class.extract_atom("mass", 2)
        logging.info("test 2 %s %s", atom_type_mass[1:10], self.me)

        # Gather the atom types, 1D int of n atoms length.
        atom_kind = lammps_class.gather_atoms("type", 0, 1)
        logging.info("test 3 %s %s", atom_kind[0:9], self.me)

        # Atom mass is indexed from 1 in lammps for some reason.
        # Create a new list rounded to the nearest mass integer
        atom_mass_type = [round(x) for x in atom_type_mass[0:ntypes + 1]]
        logging.info("test 4 %s %s", atom_mass_type[1:10], self.me)

        # Convert to masses
        final_masses = [atom_mass_type[particle] for particle in atom_kind]
        final_masses = np.array(final_masses)
        # Convert to elements
        final_elements = [ELEMENT_INDEX_MASS.get(mass, 1) for mass in final_masses]
        # logging.info("final elements test %s %s", self.me, final_elements[1:10])
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
        frame_data.arrays[PARTICLE_POSITIONS] = positions

    @property
    def interactions(self) -> List[ParticleInteraction]:
        """
        Returns a shallow copy of the current interactions.
        This is copied from the ASE example, but reduces the abstracted degree
        """
        return self.imd_server.service.active_interactions

    def add_interaction_to_ctype(self, interaction_forces: np.array, lammps_forces):
        """
        :param interaction_forces: External (user) forces
        :param lammps_forces: Internal lammps forces
        :return: Combined c_type forces
        """

        # Convert units from narupa standard of  kj/mol/nm to whatever units LAMMPS is using
        # For real units types LAMMPS uses  Kcal/mole-Angstrom 4.14
        # for kj-> Kcal and 10x for nm -> Angstrom
        interaction_forces = interaction_forces / self.force_factor
        # Flatten array into the ctype
        scatterable_array = interaction_forces.flatten()
        for idx, force in enumerate(scatterable_array):
            lammps_forces[idx] += force

    @_try_or_except
    def return_array_to_lammps(self, matrix_type: str, scatter_array, lammps_class):
        """
        Routine to return arrays to lammps
        :param matrix_type: Label for the matrix (eg. X, F, V.)
        :param scatter_array: The array to be MPI scattered
        :param lammp_class: the LAMMPS object
        """
        lammps_class.scatter_atoms(matrix_type, 1, 3, scatter_array)

    def find_unit_type(self, lammps_class):
        """
        Check the unit type collected from LAMMPS against the plank_values list and fid its index
        :return: The replaced units from the list.
        """
        plank_value = lammps_class.extract_global("hplanck", 1)
        logging.info("Plank value from lammps_internal %s ", plank_value)
        plank_type = min(range(len(PLANK_VALUES)), key=lambda i: abs(PLANK_VALUES[i] - plank_value))
        logging.info("Key detected %s", plank_type)
        return plank_type

    @_try_or_except
    def manipulate_lammps_internal_matrix(self, lammps_class, positions_3n, matrix_type):
        """
        This groups together the routines needed to return forces to lammps, is has been made general
        in case one day we and to do velocity renormalisation or another type of manipulation.

        :param lammps_class: LAMMPS class that contains all the needed routines
        :param positions_3n: Positon matrix needed to calcualte_imd_forces
        :param matrix_type: The matrix to eb scattered, usually f (forces), but could also be V (velocities)
        :return:
        """
        # Collect matrix from LAMMPS
        forces = self.manipulate_lammps_array(matrix_type, lammps_class)

        # Collect interaction vector from client on process 0
        if self.me == 0:
            interactions = self.interactions
            # Distribute to all other processors
            for i in range(1, self.nprocs):
                self.comm.send(interactions, dest=i, tag=9)

        # Collect interactions on all other processes from process 0
        if self.me > 0:
            interactions = self.comm.recv(source=0, tag=9)

        if matrix_type == 'f':
            # Create numpy arrays with the forces to be added
            energy_kjmol, forces_kjmol = calculate_imd_force(positions_3n, self.masses, interactions.values())

        self.add_interaction_to_ctype(forces_kjmol, forces)
        self.return_array_to_lammps(matrix_type, forces, lammps_class)

    def lammps_hook(self, lmp=None, comm=None):
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
                lammps_class = DummyLammps(self.n_atoms_in_dummy)
            except Exception as err:
                # Many possible reasons for LAMMPS failures so for the moment catch all
                raise Exception("Failed to load DummyLammps", err)
        else:
            # Make sure LAMMPS object is callable
            try:
                # logging.info(comm)
                lammps_class = lammps(ptr=lmp)  # , comm=comm)
            except Exception as err:
                # Many possible reasons for LAMMPS failures so for the moment catch all
                raise Exception("Failed to load LAMMPS wrapper", err)

        if self.topology_loop is True:
            units = self.find_unit_type(lammps_class)
            n_atoms = lammps_class.get_natoms()

            print("N_atoms is", n_atoms, self.me)
            units_type = LAMMPS_UNITS_CHECK.get(units, None)[0]
            distance_factor = LAMMPS_UNITS_CHECK.get(units, None)[1]
            force_factor = LAMMPS_UNITS_CHECK.get(units, None)[2]
            logging.info("units : %s %s %s %s", self.me, units_type, force_factor, distance_factor)
            self.n_atoms = n_atoms
            self.distance_factor = distance_factor
            self.units_type = units_type
            self.force_factor = force_factor
        else:
            n_atoms = self.n_atoms
            distance_factor = self.distance_factor
            units_type = self.units_type
            force_factor = self.force_factor

        self.comm.barrier()
        if self.topology_loop is True:
            atom_type, masses = self.gather_lammps_particle_types(lammps_class)
            self.masses = masses
            self.atom_type = atom_type

        # Extract the masses of the types, 1D float of home many
        # mass types were defined in input. Indexed from 1 not zero in lammps

        # Extract the position matrix
        positions = self.manipulate_lammps_array('x', lammps_class)
        # Copy the ctype array to numpy for processing
        positions_3n = np.fromiter(positions, dtype=np.float64, count=len(positions)).reshape(self.n_atoms, 3)
        positions_3n *= distance_factor

        # Collect client data and return to lammps internal arrays
        self.manipulate_lammps_internal_matrix(lammps_class, positions_3n, 'f')

        if self.me == 0:
            # Convert positions
            self.lammps_positions_to_frame_data(self.frame_data, positions)
            # Convert elements from list to frame data
            self.frame_data.arrays[PARTICLE_ELEMENTS] = self.atom_type
            # Send frame data
            self.frame_server.send_frame(self.frame_index, self.frame_data)
            self.frame_index += 1

        # Print every 100 cycles if python interpreter is still running
        # This helps ensure that everything in lammps is continuing to run
        if self.me == 0:
            self.frame_loop += 1
            if self.frame_loop == 100:
                self.frame_loop = 0

        self.topology_loop = False
