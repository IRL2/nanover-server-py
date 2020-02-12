# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.
"""
LAMMPS python integration with Narupa
This program can be run as a standalone using mock data or from within LAMMPS
using the python_invoke/fix command as demonstrated in the example LAMMPS inputs.
"""
import functools
import logging
from typing import List, Optional
import numpy as np

try:
    from lammps import lammps
except ImportError:
    logging.info('lammps failed to import', exc_info=True)

from narupa.app import NarupaImdApplication

from narupa.trajectory import FrameData
from narupa.trajectory.frame_data import PARTICLE_POSITIONS

# IMD related imports
from narupa.imd.imd_force import calculate_imd_force
from narupa.imd.particle_interaction import ParticleInteraction
from narupa.lammps.mock import MockLammps
from narupa.lammps.conversions import ELEMENT_INDEX_MASS, LAMMPS_UNITS_CHECK, PLANK_VALUES


class LammpsImd:
    """
    A class that can communicate with the LAMMPS program through
    its python interpreter. Upon initialisation, MPI is set up along with the Narupa server.
    The LAMMPS data is collected across all processors using GATHER and SCATTER routines
    that require mpi4py to respect the internal processor rank of LAMMPS.

    The variables that can currently be accessed are
    x : positions
    v : velocities
    f : forces

    This code executes as a LAMMPS fix allowing python to be executed after the forces are calculated.
    The particle positions can be extracted and the forces modified in the lammps ctypes.

    The main lammps_hook routine will check if it is being run from within LAMMPS or as a
    stand alone program and determine if it should use mock variables (manipulate_dummy_arrays)
    or ones available from within LAMMPS (manipulate_lammps_arrays).
    """
    need_to_collect_topology = True

    def __init__(self, port: int = None, address: str = "[::]"):
        """
        Items that should be initialised on instantiation of lammpsHook class
        The MPI routines are essential to stop thread issues that cause internal
        LAMMPS crashes
        """
        logging.basicConfig(level=logging.INFO)
        try:
            from mpi4py import MPI
            self.comm = MPI.COMM_WORLD
            me = MPI.COMM_WORLD.Get_rank()
            nprocs = MPI.COMM_WORLD.Get_size()
            self.me = me
            self.nprocs = nprocs

            if me == 0:
                logging.debug("MPI rank %s", me)
                logging.debug("MPI n processors %s ", nprocs)
        except ImportError as err:
            logging.error("Didn't find mpi4py %s", err)
            raise Exception("Failed to load load mpi4py, please make sure it is installed", err)

        # Start frame server, must come before MPI
        if me == 0:
            # TODO raise exception when this fails, i.e if port is blocked
            self.server_app = NarupaImdApplication.basic_server(name='LAMMPS-server', address=address, port=port)
            self.frame_service = self.server_app.frame_publisher
            self.imd_service = self.server_app.imd
            self.frame_index = 0
            self.frame_loop = 0

            try:
                self.frame_data = FrameData()
            except Exception as err:
                raise Exception("Failed to load FrameData", err)

            logging.info("Narupa Lammpshook initialised")
            logging.info("Serving on %s ", port)

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
            self.md_log_frequency = 1000

    def lammps_hook(self, lmp=None, comm=None):
        """
        lammps_hook is the main routine that is run within LAMMPS MD
        steps. It checks that the LAMMPS python wrapper is callable
        and then attempts to extract a 3N matrix of atomic data

        :param lmp: LAMMPS object data, only populated when running from within LAMMPS
        """

        if self.need_to_collect_topology is True:
            # Determine if a true lammps object is available or fall back to the test routine
            self.lammps_class = self._determine_lmp_status(comm, lmp)

            # Collect unchanging variables in the first MD step only (the topology loop)
            n_atoms, distance_factor, units_type, force_factor = self._extract_fundamental_factors(self.lammps_class)

            # Collect the lammps atom types and masses for a crude topology
            self.atom_type = self._gather_lammps_particle_types(self.lammps_class)
        else:
            # Useful to extract these in the event of MPI issues.
            n_atoms = self.n_atoms
            distance_factor = self.distance_factor
            units_type = self.units_type
            force_factor = self.force_factor

        # Extract the position matrix from lammps
        positions = self._extract_positions(self.lammps_class)

        # Collect client data and return to lammps internal arrays
        self._manipulate_lammps_internal_matrix(self.lammps_class, positions, distance_factor, 'f')

        if self.me == 0:
            # Send frame data on proc 0 only as otherwise port blocking occurs
            self._collect_and_send_frame_data(positions)

        # Print regularly if python interpreter is still running
        # This helps ensure that everything in lammps is continuing to run
        # and that the python interpreter has not crashed
        if self.me == 0:
            self.frame_loop += 1
            if self.frame_loop == self.md_log_frequency:
                self.frame_loop = 0
                logging.info("Narupa enabled calculation is still running")

        self.need_to_collect_topology = False

    @property
    def interactions(self) -> List[ParticleInteraction]:
        """
        Returns a shallow copy of the current interactions.
        This is copied from the ASE example, but reduces the abstracted degree
        """
        return self.imd_service.active_interactions

    def close(self):
        """
        Close ports to prevent blocking
        """
        logging.debug("Closing Narupa server")
        self.server_app.close()

    def _try_or_except(func):
        """
        Function creates an except or try for various functions to overcome the issue
        of the LAMMPS interpreter not giving error messages correctly when an error is
        encountered.
        :param func: function to be decorated with a try or except statement
        :return: the original function but decorated
        """

        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            result = None
            try:
                result = func(*args, **kwargs)
            except Exception as e:
                # Note args[0] is used to get around the issue of passing self to a decorator
                logging.info("Exception raised in calling function on proc %s ", args[0].me)
                logging.info("Exception thrown %s ", func)
                logging.info("Exception thrown %s ", e)
                raise
            return result

        return wrapper


    @_try_or_except
    def _gather_lammps_array(self, matrix_type: str, lammps_class):
        """
        Gather Matrix data from all LAMMPS MPI threads

        :param matrix_type: String identifying data to transmit, e.g x, v or f
        :param lammps_class: LAMMPS class that contains all the needed routines
            type.
        :return: 3N matrix with all the data requested
        """

        data_array = lammps_class.gather_atoms(matrix_type, 1, 3)

        return data_array

    def _gather_lammps_particle_types(self, lammps_class):
        """
        Collect the particle list from LAMMPS, this may be atomistic or coarse grained
        particles by looking up the ID of the atoms and that ID's corresponding mass
        from the atoms_elements dict.

        :param lammps_class: LAMMPS class that contains all the needed routines
        :return: 1N matrix with all the data requested
        """
        # Extract the number of atoms types in the system.
        ntypes = lammps_class.extract_global("ntypes", 0)

        # Extract the masses of the types, 1D float of home many
        # mass types were defined in input. Indexed from 1 not zero in lammps
        atom_type_mass = lammps_class.extract_atom("mass", 2)

        # Gather the atom types, 1D int of n atoms length.
        atom_kind = lammps_class.gather_atoms("type", 0, 1)

        # Atom mass is indexed from 1 in lammps for some reason.
        # Create a new list rounded to the nearest mass integer
        if self.units_type is "lj":
            # Some lennard jones calculations don't allow iteration over atom_type_mass
            # This is a workaround
            atom_mass_type = [0]
            for i in range(ntypes):
                atom_mass_type.append(12 + i)
        else:
            atom_mass_type = [round(x) for x in atom_type_mass[0:ntypes + 1]]

        self._log_mpi("atom_mass_types %s", atom_mass_type)
        # Convert to masses
        final_masses = [atom_mass_type[particle] for particle in atom_kind]
        final_masses = np.array(final_masses)
        # Convert to elements
        final_elements = [ELEMENT_INDEX_MASS.get(mass, 1) for mass in final_masses]
        return final_elements

    def _lammps_positions_to_frame_data(self,
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
        positions = np.ctypeslib.as_array(data_array, shape=(len(data_array),))
        # Convert to nm
        positions = positions / self.distance_factor
        self._add_pos_to_framedata(frame_data, positions)

    def _add_pos_to_framedata(self, frame_data, positions):
        frame_data.arrays[PARTICLE_POSITIONS] = positions

    def _add_interaction_to_ctype(self, interaction_forces: np.array, lammps_forces):
        """
        Adds the interaction forces to the LAMMPS array

        :param interaction_forces: External (user) forces
        :param lammps_forces: Internal lammps forces
        :return: Combined c_type forces
        """

        # Convert units from narupa standard of  kj/mol/nm to whatever units LAMMPS is using
        # For real units types LAMMPS uses  Kcal/mole-Angstrom 4.14
        # for kj-> Kcal and 10x for nm -> Angstrom
        # Flatten array into the ctype
        interaction_forces = interaction_forces.ravel()
        interaction_forces = interaction_forces / self.force_factor
        buffer = np.frombuffer(lammps_forces)
        buffer += interaction_forces

    @_try_or_except
    def _return_array_to_lammps(self, matrix_type: str, scatter_array, lammps_class):
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
        self._log_mpi("Plank value from lammps_internal %s ", plank_value)
        plank_type = min(range(len(PLANK_VALUES)), key=lambda i: abs(PLANK_VALUES[i] - plank_value))
        self._log_mpi("Key detected %s", plank_type)
        return plank_type

    @_try_or_except
    def _manipulate_lammps_internal_matrix(self, lammps_class, positions, distance_factor, matrix_type):
        """
        This groups together the routines needed to return forces to lammps,
        is has been made general in case one day we and to do velocity renormalisation
         or another type of manipulation.

        :param lammps_class: LAMMPS class that contains all the needed routines
        :param positions_3n: Position matrix needed to calculate_imd_forces
        :param matrix_type: The matrix to be scattered, usually f (forces),
            but could also be V (velocities)
        :return:
        """
        # Collect matrix from LAMMPS
        forces = self._gather_lammps_array(matrix_type, lammps_class)

        # Convert the positions to a 2D, 3N array for use in calculate)imd_force
        positions_3n = np.ctypeslib.as_array(positions, shape=(self.n_atoms * 3)).reshape(self.n_atoms, 3)

        # Collect interaction vector from client on process 0
        if self.me == 0:
            interactions = self.interactions
            # Distribute to all other processors
            for i in range(1, self.nprocs):
                self.comm.send(interactions, dest=i, tag=9)

        # Collect interactions on all other processes from process 0
        else:
            interactions = self.comm.recv(source=0, tag=9)

        if matrix_type == 'f':
            # Create numpy arrays with the forces to be added
            energy_kjmol, forces_kjmol = calculate_imd_force(positions_3n, self.masses, interactions.values())

        self._add_interaction_to_ctype(forces_kjmol, forces)
        self._return_array_to_lammps(matrix_type, forces, lammps_class)

    def _extract_positions(self, lammps_class):
        """
        Extract the particle positions and add them to the framedata.

        Return: the positions array for use in the frame data
        """

        positions = self._gather_lammps_array('x', lammps_class)

        return positions

    def _determine_lmp_status(self, comm, lmp):
        """
        Checks if the real class is being passed
        If not assume we are in interactive mode and use the DummyLammps

        return: the lammps_class object
        """
        if lmp is None:
            print("Running without lammps, assuming interactive debugging")
            try:
                lammps_class = MockLammps(self.n_atoms_in_dummy)
            except Exception as err:
                # Many possible reasons for LAMMPS failures so for the moment catch all
                raise Exception("Failed to load DummyLammps", err)
        else:
            # Make sure LAMMPS object is callable
            try:
                lammps_class = lammps(ptr=lmp, comm=comm)
            except Exception as err:
                # Many possible reasons for LAMMPS failures so for the moment catch all
                raise Exception("Failed to load LAMMPS wrapper", err)
        return lammps_class

    def _extract_fundamental_factors(self, lammps_class):
        """
        Extract the fundamental unit parameters from lammps and determine the unit types from the result
        this is only called in the first loop as the topology shouldn't change in MD simualtions after the
        initial setup.

        :param lammps_class:
        """
        units = self.find_unit_type(lammps_class)
        n_atoms = lammps_class.get_natoms()
        units_type = LAMMPS_UNITS_CHECK.get(units, None)[0]
        distance_factor = LAMMPS_UNITS_CHECK.get(units, None)[1]
        force_factor = LAMMPS_UNITS_CHECK.get(units, None)[2]
        self._log_mpi("units : %s %s %s %s", self.me, units_type, force_factor, distance_factor)
        self.n_atoms = n_atoms
        self.distance_factor = distance_factor
        self.units_type = units_type
        self.force_factor = force_factor
        return n_atoms, distance_factor, units_type, force_factor

    def _collect_and_send_frame_data(self, positions):
        """
        Collect all the data that had been extracted from lammps and transmit it in framedata

        :param positions: the flat (1D) positions arrays
        """

        self.frame_data.particle_count = self.n_atoms
        self.frame_data.particle_elements = self.atom_type
        self._lammps_positions_to_frame_data(self.frame_data, positions)

        # Send frame data
        self.frame_service.send_frame(self.frame_index, self.frame_data)
        self.frame_index += 1

    def _log_mpi(self, passed_string: str = None, *args: Optional, **kwargs: Optional):
        """
        Wrapper function for printing on one core only

        :param passed_string: string that is to be printed
        :param args: mutable objects passed
        :param kwargs: immutable objects passed
        """
        if self.me == 0:
            logging.debug(passed_string, *args, **kwargs)
