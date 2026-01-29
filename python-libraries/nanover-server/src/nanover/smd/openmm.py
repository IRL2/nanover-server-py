import os.path
from os import PathLike
from pathlib import Path
from typing import Union, Any
from abc import abstractmethod

import numpy as np

from openmm import CustomExternalForce, CustomCentroidBondForce, OpenMMException
from openmm.app import Simulation
from openmm.unit import picosecond

from nanover.openmm import serializer

SMD_FORCE_CONSTANT_PARAMETER_NAME = "smd_k"
SMD_FORCE_EXPRESSION_ATOM_NONPERIODIC = (
    f"0.5 * {SMD_FORCE_CONSTANT_PARAMETER_NAME} * ((x-x0)^2 + (y-y0)^2 + (z-z0)^2)"
)
SMD_FORCE_EXPRESSION_ATOM_PERIODIC = f"0.5 * {SMD_FORCE_CONSTANT_PARAMETER_NAME} * periodicdistance(x, y, z, x0, y0, z0)^2"
SMD_FORCE_EXPRESSION_COM = f"0.5 * {SMD_FORCE_CONSTANT_PARAMETER_NAME} * pointdistance(x1, y1, z1, x0, y0, z0)^2"

from nanover.smd.pathsmoothing import interpolate_path

from numba import jit, prange


class OpenMMSMDSimulation:
    """
    A wrapper for performing constant velocity SMD on an OpenMM simulation.

    This base class defines much of the functionality required for an SMD
    simulation to be performed using OpenMM, and automatically returns an
    instance of the appropriate subclass, depending on the mode of interaction
    required to perform SMD (i.e. either single atom or COM).
    """

    @classmethod
    def from_simulation(
        cls,
        simulation: Simulation,
        smd_atom_indices: np.ndarray,
        smd_path: np.ndarray,
        smd_force_constant: float,
        *,
        name: str | None = None,
    ):
        """
        Construct the SMD simulation from an existing OpenMM simulation.

        :param simulation: An existing OpenMM Simulation
        :param smd_atom_indices: A NumPy array of indices of the atoms to which the SMD force
          should be applied (0-D for single atom, 1-D for COM)
        :param smd_path: A NumPy array of coordinates defining the path that the
          SMD force will take during the SMD simulation
        :param smd_force_constant: The force constant of the SMD force
        :param name: An optional name for the simulation instead of default
        """

        # Create instance of SMD simulation based on type of indices passed
        assert smd_atom_indices.size >= 1
        if smd_atom_indices.size > 1:
            sim = super(cls, OpenMMSMDSimulationCOM).__new__(OpenMMSMDSimulationCOM)
        else:
            sim = super(cls, OpenMMSMDSimulationAtom).__new__(OpenMMSMDSimulationAtom)

        sim.name = name
        sim.simulation = simulation

        # Check if simulation employs periodic boundary conditions
        sim._sim_uses_pbcs = sim.simulation.system.usesPeriodicBoundaryConditions()

        # Initialise all objects relevant to the SMD simulation
        sim._initialise_smd_simulation(smd_atom_indices, smd_path, smd_force_constant)

        # Create a checkpoint of the simulation
        sim.checkpoint = sim.simulation.context.createCheckpoint()

        return sim

    @classmethod
    def from_xml_path(
        cls,
        path: PathLike[str],
        smd_atom_indices: np.ndarray,
        smd_path: np.ndarray,
        smd_force_constant: float,
        *,
        name: str | None = None,
    ):
        """
        Construct the SMD simulation from an existing NanoVer OpenMM XML file located at a given path.

        :param path: Path of the NanoVer OpenMM XML file
        :param smd_atom_indices: The indices of the atoms to which the SMD force
          should be applied
        :param smd_path: A NumPy array of coordinates defining the path that the
          SMD force will take during the SMD simulation
        :param smd_force_constant: The force constant of the SMD force
        :param name: An optional name for the simulation instead of default
        """

        if smd_atom_indices.size > 1:
            sim = super(cls, OpenMMSMDSimulationCOM).__new__(OpenMMSMDSimulationCOM)
        elif smd_atom_indices.size == 1:
            sim = super(cls, OpenMMSMDSimulationAtom).__new__(OpenMMSMDSimulationAtom)

        if name is None:
            sim.name = Path(path).stem
        else:
            sim.name = name

        sim.xml_path = path

        # Load the simulation from the path
        with open(sim.xml_path) as infile:
            sim.simulation = serializer.deserialize_simulation(infile)

        # Check if simulation employs periodic boundary conditions
        sim._sim_uses_pbcs = sim.simulation.system.usesPeriodicBoundaryConditions()

        # Initialise all objects relevant to the SMD simulation
        sim._initialise_smd_simulation(smd_atom_indices, smd_path, smd_force_constant)

        # Create a checkpoint of the simulation
        sim.checkpoint = sim.simulation.context.createCheckpoint()

        return sim

    def __init__(self, name: str | None = None):
        self.name = name or "Unnamed OpenMM SMD Simulation"

        self.xml_path: PathLike[str] | None = None

        self.simulation: Simulation | None = None
        self.smd_atom_indices: np.ndarray | None = None
        self.smd_path: np.ndarray | None = None
        self.smd_force_constant: float | None = None

        self.loaded_smd_force_from_sim: bool = False
        self.n_smd_atom_indices: int | None = None

        self.smd_force: Union[CustomExternalForce, CustomCentroidBondForce] | None = (
            None
        )

        self.checkpoint: Any | None = None

        self.current_smd_force_position: np.ndarray | None = None
        self.current_smd_force_position_index: int | None = None
        self.smd_simulation_atom_positions: np.ndarray | None = None
        self.smd_simulation_forces: np.ndarray | None = None
        self.smd_simulation_work_done: np.ndarray | None = None

        self._sim_uses_pbcs: bool | None = None

    @abstractmethod
    def define_smd_simulation_atom_positions_array(self):
        """
        Define the array to which the positions of the atoms with which the
        SMD force interacts will be saved over the course of the SMD simulation.
        """
        pass

    @abstractmethod
    def check_for_existing_smd_force(self):
        """
        Check whether the loaded simulation already contains an SMD force of
        the correct type for the simulation type.
        """
        pass

    @abstractmethod
    def add_smd_force_to_system(self):
        """
        Add the required SMD force to the system, depending on the type of
        SMD interaction required (single atom or centre-of-mass).
        """
        pass

    @abstractmethod
    def update_smd_force_position(self):
        """
        Update the position of the SMD force.
        """
        pass

    @abstractmethod
    def calculate_cumulative_work_done(self):
        """
        Calculate the cumulative work done along the reaction coordinate
        during the SMD simulation.
        """
        pass

    def _initialise_smd_simulation(
        self,
        smd_atom_indices: np.ndarray,
        smd_path: np.ndarray,
        smd_force_constant: float,
    ):
        """
        Set the fields relevant to the SMD simulation and add the SMD force to the system.
        Called upon when constructing an OpenMMSMDSimulation using the classmethods
        from_simulation or from_xml_path.
        """
        self.smd_atom_indices = smd_atom_indices
        self.n_smd_atom_indices = self.smd_atom_indices.size
        self.smd_path = smd_path
        self.smd_force_constant = smd_force_constant
        self.define_smd_simulation_atom_positions_array()

        self.current_smd_force_position = self.smd_path[0]
        self.current_smd_force_position_index = 0

        # Check whether SMD force is already present
        self.check_for_existing_smd_force()

        if not self.loaded_smd_force_from_sim:
            # Create SMD force and add it to the system
            self.add_smd_force_to_system()

        self.simulation.context.reinitialize(preserveState=True)

    def reset(self):
        """
        Reset the SMD simulation to its initial state, and reset the arrays output
        by the SMD simulation.
        """
        assert (
            self.simulation is not None
            and self.checkpoint is not None
            and self.smd_path is not None
            and self.smd_force is not None
        )

        # Reset SMD force position
        self.current_smd_force_position = self.smd_path[0]
        self.current_smd_force_position_index = 0
        self.update_smd_force_position()

        # Reset simulation, reinitialise context to be safe
        self.simulation.context.reinitialize()
        self.simulation.context.loadCheckpoint(self.checkpoint)

        # Reset or remove SMD simulation arrays
        self.define_smd_simulation_atom_positions_array()
        del self.smd_simulation_forces
        del self.smd_simulation_work_done

    def remove_smd_force_from_system(self):
        """
        Remove any SMD forces from the system.
        """
        forces = self.simulation.system.getForces()
        forces_to_remove = []
        for i in range(len(forces)):
            if (type(forces[i]) == type(self.smd_force)):
                try:
                    if forces[i].getGlobalParameterName(0) == SMD_FORCE_CONSTANT_PARAMETER_NAME:
                        forces_to_remove.append(i)
                except OpenMMException:
                    continue
                # forces_to_remove.append(i)

        # Remove any SMD forces, accounting for the changes in indices
        # as forces are removed
        for j in range(len(forces_to_remove)):
            self.simulation.system.removeForce(forces_to_remove[j] - j)

    def get_smd_atom_positions(self):
        """
        Retrieve the positions of the atoms with which the SMD force is
        interacting, and add them to the array of positions to save.
        """
        positions = self.simulation.context.getState(getPositions=True).getPositions(
            asNumpy=True
        )
        self.smd_simulation_atom_positions[self.current_smd_force_position_index] = (
            positions[self.smd_atom_indices]
        )

    def run_equilibration_with_initial_restraint(self, n_steps: int):
        """
        Perform an equilibration of the system with the restraint fixed at the
        initial position defined in the SMD force path.

        :param n_steps: Number of simulation steps to run equilibration for.
        """

        try:
            assert np.all(self.current_smd_force_position == self.smd_path[0])
        except AssertionError:
            raise AssertionError(
                "Restraint is not located at the initial position "
                "of the SMD force, equilibration aborted. To prepare the system "
                "appropriately, the restraint should be placed at the first point "
                "on the path that the SMD force will take during the SMD simulation."
            )

        self.simulation.step(n_steps)

    def save_simulation(
        self,
        output_filepath: PathLike[str],
        save_state: bool = False,
        save_smd_force: bool | None = False,
    ):
        """
        Save the simulation to a NanoVer OpenMM XML file, with the option to include the
        SMD force in the XML file.

        :param output_filepath: Path to output file to save the simulation to.
        :param save_state: If True, save the present state of the simulation to the XML file.
        :param save_smd_force: Bool defining whether to save the SMD force in the XML file (Optional).
        """
        assert output_filepath is not None

        if save_smd_force:
            with open(output_filepath, "w") as outfile:
                outfile.write(
                    serializer.serialize_simulation(
                        self.simulation, save_state=save_state
                    )
                )

        else:
            # Temporarily remove SMD force from simulation
            self.remove_smd_force_from_system()

            # Save simulation without SMD force
            with open(output_filepath, "w") as outfile:
                xml_string = serializer.serialize_simulation(
                    self.simulation, save_state=save_state
                )
                # Manually remove parameters for now...
                xml_string_lines = xml_string.splitlines()
                for line in xml_string_lines:
                    # Remove only parameter left over from SMD force
                    if SMD_FORCE_CONSTANT_PARAMETER_NAME in line:
                        line_index = xml_string_lines.index(line)
                        n_tabs = line.index("<")
                        param_string = line.split('/')
                        param_string_list = sum([string.split() for string in param_string], [])
                        for substring in param_string_list:
                            if SMD_FORCE_CONSTANT_PARAMETER_NAME in substring:
                                param_string_list.remove(substring)
                        xml_string_lines[line_index] = n_tabs * "\t" + " ".join(param_string_list[:-1]) + "/" + \
                                                       param_string_list[-1]
                xml_string = "\n".join(xml_string_lines)
                outfile.write(xml_string)

            # Add SMD force back to the system
            self.add_smd_force_to_system()

    def generate_starting_structures(
        self,
        interval_ps: float,
        n_structures: int,
        output_directory: PathLike[str] | None = None,
        filename_prefix: str | None = None,
        save_smd_force: bool | None = None,
    ):
        """
        Generate the specified number of starting structures by running the simulation for the specified
        interval with the initial restraint applied to the system and saving structures at regular intervals.
        Structures are saved to the output path, optionally with the filename prefix. If no output path is
        specified, structures will be saved to the current working directory. If no prefix is specified,
        the output files will be named automatically.

        :param interval_ps: Interval to run the simulation for, in picoseconds.
        :param n_structures: Number of structures to generate.
        :param output_directory: Output directory to save the structures to (Optional).
        :param filename_prefix: Prefix for output files (Optional).
        :param save_smd_force: Bool defining whether to save the SMD force in the XML file (Optional).
        """
        timestep_ps = self.simulation.integrator.getStepSize()._value
        n_steps_struct_interval = int(
            np.floor(interval_ps / (n_structures * timestep_ps))
        )

        if output_directory is None:
            output_directory = os.getcwd()

        if filename_prefix is None:
            filename_prefix = "smd_structure"

        print(
            f"Generating {n_structures} structures in {interval_ps} ps simulation...\n"
            f"Structures will be saved to {output_directory}\n"
        )

        for i in range(n_structures):

            # Run set of simulation steps to generate next structure
            self.simulation.step(n_steps_struct_interval)

            # Save structure to output file
            outfile_path = os.path.join(
                output_directory, filename_prefix + "_" + str(i + 1) + ".xml"
            )
            self.save_simulation(
                output_filepath=outfile_path,
                save_state=True,
                save_smd_force=save_smd_force,
            )

        print(f"Structure generation complete: {n_structures} structures generated.")

    def run_smd(self, progress_interval: int | None = 100000):
        """
        Perform an SMD simulation on the system using the SMD force and
        path defined, store the positions of the atoms with which the
        SMD force interacts, and calculate the work done along the reaction
        coordinate.

        :param progress_interval: Interval defining how regularly to print the progress
          of the simulation, in terms of the number of simulation steps.
        """

        # Get initial atom positions for SMD force at initial position
        assert (
            self.smd_simulation_atom_positions is not None
            and self.smd_force is not None
            and np.all(self.current_smd_force_position == self.smd_path[0])
            and self.current_smd_force_position_index == 0
        )
        self.get_smd_atom_positions()

        # Run SMD procedure
        n_steps = self.smd_path.shape[0]

        print("Starting SMD simulation...")
        for step in range(1, n_steps):

            # Perform single simulation step
            self.simulation.step(1)

            # Update force position index
            self.current_smd_force_position_index = step

            # Update SMD force position
            self.update_smd_force_position()

            # Retrieve atom positions on step
            self.get_smd_atom_positions()

            # Print step at intervals
            if step % progress_interval == 0:
                print(f"Steps completed: {step}")

        print("SMD simulation completed. Calculating work done...")

        # Calculate the work done along the reaction coordinate
        self.calculate_cumulative_work_done()

        print("Work done calculated.")

    def _calculate_smd_forces(self, interaction_centre_positions):
        """
        Calculate the SMD forces that acted on the system during the simulation in kJ mol-1 nm-1.

        :param interaction_centre_positions: Array of positions defining the centre
          (single atom or COM of group of atoms) with which the SMD force interacted during the simulation.
        """
        assert np.all(self.smd_path.shape == interaction_centre_positions.shape)
        self.smd_simulation_forces = -self.smd_force_constant * (
            interaction_centre_positions - self.smd_path
        )

    def _calculate_work_done(self):
        """
        Calculate the work done along the reaction coordinate in kJ mol-1.
        """
        assert self.smd_path is not None and self.smd_simulation_forces is not None
        smd_force_displacements = np.diff(self.smd_path, axis=0)
        work_done_array = np.zeros(self.smd_simulation_forces.shape[0])
        for i in range(smd_force_displacements.shape[0]):
            work_done_array[i + 1] = np.dot(
                self.smd_simulation_forces[i], smd_force_displacements[i]
            )
        self.smd_simulation_work_done = np.cumsum(work_done_array, axis=0)

    def save_smd_simulation_data(self,
                                 path: PathLike[str] = None,
                                 save_work_done: bool = True,
                                 save_atom_positions: bool = True,
                                 work_done_dtype: np.dtype = np.float32,
                                 atom_positions_dtype: np.dtype = np.float32):
        """
        Save the data produced by the SMD simulation in binary form that can be read
        into NumPy arrays. The following data are saved, in the order listed below:

        - Trajectories of the atoms to which the SMD force was applied, in nm
        - Work done along the reaction coordinate defined by the path of the SMD force, in kJ mol-1

        :param path: Path to the file to which the data will be saved.
        :param save_work_done: Bool determining whether to save the work done
        :param save_atom_positions: Bool determining whether to save the positions of the atom(s)
        :param work_done_dtype: Data type of the work done array to save.
        :param atom_positions_dtype: Data type of the atom positions array to save.
        """

        if path is None:
            raise ValueError(
                "Output file path cannot be None. Please specify an output file path."
            )

        elif self.smd_simulation_work_done is None:
            raise ValueError(
                "Missing values for the work done. This data can only be saved after"
                "the SMD calculation is completed."
            )

        elif self.smd_simulation_atom_positions is None or np.all(
            self.smd_simulation_atom_positions == 0.0
        ):
            raise ValueError(
                "Missing values for the atom positions. This data can only be saved after"
                "the SMD calculation is completed."
            )

        with open(path, "wb") as outfile:
            if save_atom_positions:
                np.save(outfile, self.smd_simulation_atom_positions.astype(atom_positions_dtype))
                print("Atom positions saved to simulation data file.")
            if save_work_done:
                np.save(outfile, self.smd_simulation_work_done.astype(work_done_dtype))
                print("Work done saved to simulation data file.")

    def save_general_smd_data(self, path: PathLike[str] = None):
        """
        Saves general data related to the SMD simulation in binary form that can be read
        into NumPy arrays. The following data are saved, in the order listed below:

        - Indices of the atoms to which the SMD force is applied
        - Positions defining the path of the SMD force, in nm
        - Force constant of the SMD force, in kJ mol-1 nm-2
        - Temperature of the simulation, in Kelvin
        - Time step of the simulation, in picoseconds

        :param path: Path to the file to which the data will be saved.
        """

        if path is None:
            raise ValueError(
                "Output file path cannot be None. Please specify an output file path."
            )

        with open(path, "wb") as outfile:
            np.save(outfile, self.smd_atom_indices)
            np.save(outfile, self.smd_path)
            np.save(outfile, self.smd_force_constant)
            np.save(outfile, self.simulation.integrator.getTemperature()._value)
            np.save(outfile, self.simulation.integrator.getStepSize()._value)


class OpenMMSMDSimulationAtom(OpenMMSMDSimulation):
    """
    A class for performing constant velocity SMD on an OpenMM simulations,
    where the SMD force is applied to a single atom.
    """

    def __init__(self, name: str | None = None):

        super().__init__(name)

    def check_for_existing_smd_force(self):

        try:
            force_constant = self.simulation.context.getParameter(
                SMD_FORCE_CONSTANT_PARAMETER_NAME
            )
            n_forces = self.simulation.system.getNumForces()
            smd_force = self.simulation.system.getForce(n_forces - 1)
            params = smd_force.getParticleParameters(0)
            assert type(smd_force) == CustomExternalForce
            assert (
                smd_force.getEnergyFunction() == SMD_FORCE_EXPRESSION_ATOM_PERIODIC
                or smd_force.getEnergyFunction()
                == SMD_FORCE_EXPRESSION_ATOM_NONPERIODIC
            )
            assert smd_force.getNumParticles() == 1
            assert force_constant == self.smd_force_constant
            assert params[0] == self.smd_atom_indices
            assert np.all(params[1] == self.current_smd_force_position)
            print("SMD force already present in loaded simulation.")
            self.smd_force = smd_force
            self.loaded_smd_force_from_sim = True

        except OpenMMException:
            self.loaded_smd_force_from_sim = False

    def add_smd_force_to_system(self):
        """
        Add an SMD force to the OpenMM system that interacts with the
        specified atom.
        """

        x0, y0, z0 = self.current_smd_force_position
        smd_force = smd_single_atom_force(self.smd_force_constant, self._sim_uses_pbcs)
        smd_force.addParticle(self.smd_atom_indices, [x0, y0, z0])
        self.smd_force = smd_force
        self.simulation.system.addForce(self.smd_force)
        self.simulation.context.reinitialize(preserveState=True)

    def update_smd_force_position(self):

        self.current_smd_force_position = self.smd_path[
            self.current_smd_force_position_index
        ]
        x0, y0, z0 = self.current_smd_force_position
        self.smd_force.setParticleParameters(0, self.smd_atom_indices, [x0, y0, z0])
        self.smd_force.updateParametersInContext(self.simulation.context)

    def define_smd_simulation_atom_positions_array(self):
        self.smd_simulation_atom_positions = np.zeros((self.smd_path.shape[0], 3))

    def calculate_cumulative_work_done(self):
        """
        Calculate the cumulative work done by the SMD force on the
        atom with which it interacts over the SMD simulation.
        """
        assert not np.array_equal(
            self.smd_simulation_atom_positions, np.zeros((self.smd_path.shape[0], 3))
        )
        self._calculate_smd_forces(self.smd_simulation_atom_positions)
        self._calculate_work_done()


class OpenMMSMDSimulationCOM(OpenMMSMDSimulation):
    """
    A class for performing constant velocity SMD on an OpenMM simulations,
    where the SMD force is applied to the centre of mass of a specified
    group of atoms.
    """

    def __init__(self, name: str | None = None):

        super().__init__(name)
        self.com_positions: np.ndarray | None = None

    def check_for_existing_smd_force(self):
        try:
            force_constant = self.simulation.context.getParameter(
                SMD_FORCE_CONSTANT_PARAMETER_NAME
            )
            n_forces = self.simulation.system.getNumForces()
            smd_force = self.simulation.system.getForce(n_forces - 1)
            assert type(smd_force) == CustomCentroidBondForce
            assert smd_force.getEnergyFunction() == SMD_FORCE_EXPRESSION_COM
            assert smd_force.getNumGroups() == 1
            assert force_constant == self.smd_force_constant
            assert np.all(
                (smd_force.getGroupParameters(0)[0] == self.smd_atom_indices) == True
            )
            assert np.all(
                smd_force.getBondParameters(0)[1] == self.current_smd_force_position
            )
            print("SMD force already present in loaded simulation.")
            self.smd_force = smd_force
            self.loaded_smd_force_from_sim = True

        except OpenMMException:
            self.loaded_smd_force_from_sim = False

    def add_smd_force_to_system(self):
        """
        Add an SMD force to the OpenMM system that interacts with the
        centre of mass of the specified group of atoms.
        """

        x0, y0, z0 = self.current_smd_force_position
        smd_force = smd_com_force(self.smd_force_constant, self._sim_uses_pbcs)
        smd_force.addGroup(self.smd_atom_indices)
        smd_force.addBond([0], [x0, y0, z0])
        self.smd_force = smd_force
        self.simulation.system.addForce(self.smd_force)
        self.simulation.context.reinitialize(preserveState=True)

    def update_smd_force_position(self):

        self.current_smd_force_position = self.smd_path[
            self.current_smd_force_position_index
        ]
        x0, y0, z0 = self.current_smd_force_position
        self.smd_force.setBondParameters(0, [0], [x0, y0, z0])
        self.smd_force.updateParametersInContext(self.simulation.context)

    def define_smd_simulation_atom_positions_array(self):
        self.smd_simulation_atom_positions = np.zeros(
            (self.smd_path.shape[0], self.smd_atom_indices.size, 3)
        )

    def _calculate_com(self, atom_positions: np.ndarray, atom_masses: np.ndarray):
        """
        Calculate the centre of mass of the N atoms to which the SMD force has been applied.

        :param atom_positions: NumPy array of atom positions with dimensions (N, 3)
        :param atom_masses: NumPy array of atomic masses (AMU) with dimension (N)
        :return: NumPy array containing the position of the centre of mass of the atoms with dimension (3)
        """
        assert np.all(atom_positions.shape == np.array((self.smd_atom_indices.size, 3)))
        return calculate_com(atom_positions, atom_masses)

    def _calculate_com_trajectory(self):
        """
        Calculate the trajectory that the COM follows during the SMD simulation.
        """
        assert not np.array_equal(
            self.smd_simulation_atom_positions,
            np.zeros((self.smd_path.shape[0], self.smd_atom_indices.size, 3)),
        )
        atom_masses = np.zeros(self.n_smd_atom_indices)
        for index in range(self.n_smd_atom_indices):
            atom_masses[index] = self.simulation.system.getParticleMass(
                self.smd_atom_indices[index]
            )._value
        self.com_positions = np.array(
            [
                self._calculate_com(self.smd_simulation_atom_positions[i], atom_masses)
                for i in range(self.smd_path.shape[0])
            ]
        )

    def calculate_cumulative_work_done(self):
        """
        Calculate the cumulative work done by the SMD force on the
        COM of the atoms with which it interacts over the SMD simulation.
        """
        self._calculate_com_trajectory()
        self._calculate_smd_forces(self.com_positions)
        self._calculate_work_done()

    def save_smd_simulation_data(self,
                                 path: PathLike[str] = None,
                                 save_com_positions: bool = True,
                                 com_positions_dtype: np.dtype = np.float32,
                                 **kwargs):

        """
        Save the data produced by the SMD simulation in binary form that can be read
        into NumPy arrays. The following data can be saved (optionally), in the order listed below:

        - Trajectories of the atoms to which the SMD force was applied, in nm
        - Work done along the reaction coordinate defined by the path of the SMD force, in kJ mol-1
        - Trajectories of the COM of the atoms to which the SMD force was applied, in nm

        :param path: Path to the file to which the data will be saved.
        :param save_work_done: Bool determining whether to save the work done
        :param save_atom_positions: Bool determining whether to save the positions of the atom(s)
        :param save_com_positions: Bool determining whether to save the positions of the COM
        :param work_done_dtype: Data type of the work done array to save.
        :param atom_positions_dtype: Data type of the atom positions array to save.
        :param com_positions_dtype: Data type of the COM positions array to save.
        """
        #TODO: Think about whether there is a cleaner way to achieve this

        # Optionally save work done and atomic coordinates
        super().save_smd_simulation_data(path, **kwargs)

        if self.com_positions is None or np.all(
            self.com_positions == 0.0
        ):
            raise ValueError(
                "Missing values for the atom positions. This data can only be saved after"
                "the SMD calculation is completed."
            )

        # Optionally save COM positions
        with open(path, "ab+") as outfile:
            if save_com_positions:
                np.save(outfile, self.com_positions.astype(com_positions_dtype))
                print("COM positions saved to simulation data file.")


class OpenMMStringOptimiser:

    """
    Class for optimising the path defining the trajectory of an SMD force
    during an SMD simulation using the finite temperature string method.
    """

    def __init__(self, smd_simulation: Union[OpenMMSMDSimulationAtom, OpenMMSMDSimulationCOM], **kwargs):

        self.smd_simulation = smd_simulation
        self.checkpoints: Optional[Any] = None
        self.node_positions: Optional[np.ndarray] = None
        self.tangent_vectors: Optional[np.ndarray] = None

        self.node_position_history = []

    def generate_nodes(self, n_nodes):
        """
        Generate the checkpoint files defining the starting node positions and
        contexts for the string optimisation.

        :param n_nodes: Number of nodes to generate along the SMD force path (excluding
          the zeroth node marking the initial geometry).
        """

        # Calculate number of steps to run between nodes, discarding remainder
        # to produce n_nodes evenly spaced nodes (except ends)
        n_smd_steps = self.smd_simulation.smd_path.shape[0]
        node_interval = n_smd_steps//n_nodes

        checkpoints = []
        node_positions = np.zeros((n_nodes+1, 3))

        # Save initial checkpoint and force position
        assert np.all(self.smd_simulation.current_smd_force_position == self.smd_simulation.smd_path[0])
        checkpoints.append(self.smd_simulation.simulation.context.createCheckpoint())
        node_positions[0, :] = self.smd_simulation.current_smd_force_position

        print("Generated checkpoint 0 (initial force position)")

        for i in range(n_nodes):

            step_base = i * node_interval

            for step in range(node_interval):

                force_position_index = step + step_base
                #print(f"Step {force_position_index}")

                # Update force position index
                self.smd_simulation.current_smd_force_position_index = force_position_index

                #print(self.smd_simulation.smd_path[self.smd_simulation.current_smd_force_position_index])

                # Update SMD force position
                self.smd_simulation.update_smd_force_position()

                #print(self.smd_simulation.current_smd_force_position)

                # Perform single simulation step
                self.smd_simulation.simulation.step(1)

            # Save checkpoint and force position
            checkpoints.append(self.smd_simulation.simulation.context.createCheckpoint())
            node_positions[i+1, :] = self.smd_simulation.current_smd_force_position

            # print(f"Steps completed: {force_position_index}")
            print(f"Generated checkpoint {i+1}")

        print("All checkpoints generated")

        self.checkpoints = checkpoints
        self.node_positions = node_positions
        self.node_position_history.append(self.node_positions)
        # self.smd_simulation.reset()

    def calculate_tangent_vectors(self):

        assert self.node_positions is not None
        n_nodes = self.node_positions.shape[0]

        tangent_vectors = np.zeros((n_nodes, 3))

        # Determine tangent vectors for all but start and end nodes
        differences = np.diff(self.node_positions, axis=0)
        for i in range(1, n_nodes - 1):
            tangent_vec = differences[i-1] + differences[i]
            normalised_tangent_vec = tangent_vec / np.linalg.norm(tangent_vec)
            tangent_vectors[i] = normalised_tangent_vec

        # Determine tangent vectors for start and end nodes
        tangent_vectors[0] = differences[0] / np.linalg.norm(differences[0])
        tangent_vectors[-1] = differences[-1] / np.linalg.norm(differences[-1])

        self.tangent_vectors = tangent_vectors


    def calculate_com(self, atom_positions: np.ndarray, atom_masses: np.ndarray):
        """
        Calculate the centre of mass of a group of N atoms.

        :param atom_positions: NumPy array of atom positions with dimensions (N, 3)
        :param atom_masses: NumPy array of atomic masses (AMU) with dimension (N)
        :return: NumPy array containing the position of the centre of mass of the atoms with dimension (3)
        """
        # TODO: Check this works and write tests!
        assert np.all(atom_positions.shape == np.array((self.smd_simulation.smd_atom_indices.size, 3)))
        return np.sum(
            np.multiply(np.transpose(atom_positions), atom_masses), axis=1
        ) / np.sum(atom_masses)

    def update_force_position(self, position: np.ndarray):

        x0, y0, z0 = position
        self.smd_simulation.smd_force.setParticleParameters(0, self.smd_simulation.smd_atom_indices, [x0, y0, z0])
        self.smd_simulation.smd_force.updateParametersInContext(self.smd_simulation.simulation.context)

    def components_of_vector_parallel_to_RC(self, distance: np.ndarray, tangent: np.ndarray):
        inner_force_velocity = np.dot(distance, tangent)
        norm_squared_velocity = np.dot(tangent, tangent)
        parallel_component = (inner_force_velocity / norm_squared_velocity) * tangent
        return parallel_component

    def components_of_vector_perpendicular_to_RC_old(self, distance: np.ndarray, tangent: np.ndarray):
        parallel_component = self.components_of_vector_parallel_to_RC(distance, tangent)
        perpendicular_component = distance - parallel_component
        return perpendicular_component

    def components_of_vector_perpendicular_to_RC(self, vector: np.ndarray, tangent: np.ndarray):
        """
        Calculate the perpendicular component of a vector with respect to the reaction coordinate.
        """
        return vector - np.dot(vector, tangent) * tangent

    def update_timestep_ps(self,  timestep_ps: float):
        for i in range(len(self.checkpoints)):
            self.smd_simulation.simulation.context.loadCheckpoint(self.checkpoints[i])
            self.smd_simulation.simulation.integrator.setStepSize(timestep_ps*picosecond)
            self.smd_simulation.simulation.context.reinitialize()
            self.checkpoints[i] = self.smd_simulation.simulation.context.createCheckpoint()

        # TODO: Write test to confirm this works as expected

    def equilibrate_checkpoints(self, n_steps: int):

        assert self.checkpoints is not None

        for i in range(len(self.checkpoints)):

            self.smd_simulation.simulation.context.reinitialize()
            self.smd_simulation.simulation.context.loadCheckpoint(self.checkpoints[i])

            self.smd_simulation.simulation.step(n_steps)

            self.checkpoints[i] = self.smd_simulation.simulation.context.createCheckpoint()

    def smooth_path_naively(self, smoothing_parameter: float):
        """
        Naively smooth path without changing end points

        :param smoothing_parameter: Float between 0 and 1 defining the smoothing parameter to
          use for the path smoothing
        """
        assert self.node_positions is not None
        self.node_positions[1:-1] = (1-smoothing_parameter) * self.node_positions[1:-1] + 0.5 * smoothing_parameter * (self.node_positions[0:-2]+self.node_positions[2:])

    def calculate_minimum_energy_path(self, n_iterations: int, smoothing_parameter: float, fix_initial_position: bool = True, fix_final_position: bool = True, max_update_distance_nm: float = 0.005):

        assert self.checkpoints is not None

        atom_masses = np.zeros(self.smd_simulation.n_smd_atom_indices)
        for index in range(self.smd_simulation.n_smd_atom_indices):
            atom_masses[index] = self.smd_simulation.simulation.system.getParticleMass(
                self.smd_simulation.smd_atom_indices[index]
            )._value

        if fix_initial_position and fix_final_position:
            j_range = prange(1, len(self.checkpoints) - 1)
        elif fix_initial_position and not fix_final_position:
            j_range = prange(1, len(self.checkpoints))
        elif fix_final_position and not fix_initial_position:
            j_range = prange(len(self.checkpoints) - 1)
        else:
            j_range = prange(len(self.checkpoints))


        self.calculate_tangent_vectors()

        for j in j_range:
            self.smd_simulation.simulation.context.reinitialize()
            self.smd_simulation.simulation.context.loadCheckpoint(self.checkpoints[j])
            self.update_force_position(self.node_positions[j])
            self.checkpoints[j] = self.smd_simulation.simulation.context.createCheckpoint()

        for i in range(n_iterations):

            print(f"Iteration {i+1}")

            self.calculate_tangent_vectors()

            for j in j_range:

                self.smd_simulation.simulation.context.reinitialize()
                self.smd_simulation.simulation.context.loadCheckpoint(self.checkpoints[j])

                self.smd_simulation.simulation.minimizeEnergy()

                positions = self.smd_simulation.simulation.context.getState(getPositions=True).getPositions(
                    asNumpy=True
                )
                positions_array = (
                    positions[self.smd_simulation.smd_atom_indices]
                )

                # Calculate new position of restraint from mean of COM
                com_position = self.calculate_com(positions_array, atom_masses)
                perpendicular_component = self.components_of_vector_perpendicular_to_RC(com_position - self.node_positions[j], self.tangent_vectors[j])
                # Scale the update based on some threshold, commented out for now
                perpendicular_component_magnitude = np.linalg.norm(perpendicular_component)
                update_distance = np.min([perpendicular_component_magnitude, max_update_distance_nm])
                normalised_perpendicular_vector = perpendicular_component / perpendicular_component_magnitude
                update_vector = update_distance * normalised_perpendicular_vector
                new_node_position = self.node_positions[j] + update_vector
                #new_node_position = self.node_positions[j] + perpendicular_component

                # Update position of force and save new checkpoint
                self.update_force_position(new_node_position)
                self.checkpoints[j] = self.smd_simulation.simulation.context.createCheckpoint()

                self.node_positions[j] = new_node_position


            self.smooth_path_naively(smoothing_parameter=smoothing_parameter)
            self.node_positions, _ = interpolate_path(self.node_positions[:, 0],
                                               self.node_positions[:, 1],
                                               self.node_positions[:, 2],
                                               0.0,
                                               self.node_positions.shape[0])

            self.node_position_history.append(self.node_positions)

            # Add evenly spaced node positions to checkpoints
            for j in j_range:
                self.smd_simulation.simulation.context.reinitialize()
                self.smd_simulation.simulation.context.loadCheckpoint(self.checkpoints[j])
                self.update_force_position(self.node_positions[j])
                self.checkpoints[j] = self.smd_simulation.simulation.context.createCheckpoint()

        print(f"Completed {n_iterations} iterations.")

    # @jit(parallel=True)
    def optimise_string(self, n_sim_steps: int, n_iterations: int, smoothing_parameter: float, fix_initial_position: bool = True, fix_final_position: bool = True, max_update_distance_nm: float = 0.005):

        assert self.checkpoints is not None

        atom_masses = np.zeros(self.smd_simulation.n_smd_atom_indices)
        for index in range(self.smd_simulation.n_smd_atom_indices):
            atom_masses[index] = self.smd_simulation.simulation.system.getParticleMass(
                self.smd_simulation.smd_atom_indices[index]
            )._value

        if fix_initial_position and fix_final_position:
            j_range = prange(1, len(self.checkpoints) - 1)
        elif fix_initial_position and not fix_final_position:
            j_range = prange(1, len(self.checkpoints))
        elif fix_final_position and not fix_initial_position:
            j_range = prange(len(self.checkpoints) - 1)
        else:
            j_range = prange(len(self.checkpoints))


        self.calculate_tangent_vectors()

        for j in j_range:
            self.smd_simulation.simulation.context.reinitialize()
            self.smd_simulation.simulation.context.loadCheckpoint(self.checkpoints[j])
            self.update_force_position(self.node_positions[j])
            self.checkpoints[j] = self.smd_simulation.simulation.context.createCheckpoint()

        for i in range(n_iterations):

            print(f"Iteration {i+1}")

            self.calculate_tangent_vectors()

            for j in j_range:
                print(j)

                self.smd_simulation.simulation.context.reinitialize()
                self.smd_simulation.simulation.context.loadCheckpoint(self.checkpoints[j])

                positions_array = np.zeros((n_sim_steps, self.smd_simulation.n_smd_atom_indices, 3))

                for k in range(n_sim_steps):

                    self.smd_simulation.simulation.step(1)

                    positions = self.smd_simulation.simulation.context.getState(getPositions=True).getPositions(
                        asNumpy=True
                    )
                    positions_array[k, :, :] = (
                        positions[self.smd_simulation.smd_atom_indices]
                    )

                # Calculate new position of restraint from mean of COM
                mean_COM_position = np.mean([self.calculate_com(positions_array[l], atom_masses) for l in range(n_sim_steps)], axis=0)
                perpendicular_component = self.components_of_vector_perpendicular_to_RC(mean_COM_position - self.node_positions[j], self.tangent_vectors[j])
                # Scale the update based on some threshold, commented out for now
                perpendicular_component_magnitude = np.linalg.norm(perpendicular_component)
                update_distance = np.min([perpendicular_component_magnitude, max_update_distance_nm])
                normalised_perpendicular_vector = perpendicular_component / perpendicular_component_magnitude
                update_vector = update_distance * normalised_perpendicular_vector
                new_node_position = self.node_positions[j] + update_vector
                #new_node_position = self.node_positions[j] + perpendicular_component

                # Update position of force and save new checkpoint
                self.update_force_position(new_node_position)
                self.checkpoints[j] = self.smd_simulation.simulation.context.createCheckpoint()

                self.node_positions[j] = new_node_position


            self.smooth_path_naively(smoothing_parameter=smoothing_parameter)
            self.node_positions, _ = interpolate_path(self.node_positions[:, 0],
                                               self.node_positions[:, 1],
                                               self.node_positions[:, 2],
                                               0.0,
                                               self.node_positions.shape[0])

            self.node_position_history.append(self.node_positions)

            # Add evenly spaced node positions to checkpoints
            for j in j_range:
                self.smd_simulation.simulation.context.reinitialize()
                self.smd_simulation.simulation.context.loadCheckpoint(self.checkpoints[j])
                self.update_force_position(self.node_positions[j])
                self.checkpoints[j] = self.smd_simulation.simulation.context.createCheckpoint()

        print(f"Completed {n_iterations} iterations.")



def smd_com_force(force_constant: float, uses_pbcs: bool):
    """
    Defines a harmonic restraint force for the COM of a group of atoms for performing SMD.

    :param force_constant: Force constant of the harmonic restraint to be applied to the
      COM of the group of atoms in units kJ mol-1 nm-2
    :param uses_pbcs: Bool specifying whether to use periodic boundary conditions for the
      harmonic restraint
    :return: CustomCentroidBondForce defining the harmonic SMD force that interacts with
      the COM of the specified atoms.
    """

    smd_force = CustomCentroidBondForce(1, SMD_FORCE_EXPRESSION_COM)
    smd_force.addGlobalParameter(SMD_FORCE_CONSTANT_PARAMETER_NAME, force_constant)
    smd_force.addPerBondParameter("x0")
    smd_force.addPerBondParameter("y0")
    smd_force.addPerBondParameter("z0")
    smd_force.setUsesPeriodicBoundaryConditions(uses_pbcs)
    smd_force.setForceGroup(31)

    return smd_force


def smd_single_atom_force(force_constant: float, uses_pbcs: bool):
    """
    Defines a harmonic restraint force for a single atom for performing SMD.

    :param force_constant: Force constant of the harmonic restraint to be applied to the
      specified atom in units kJ mol-1 nm-2
    :param uses_pbcs: Bool specifying whether to use periodic boundary conditions for the
      harmonic restraint
    :return: CustomExternalForce defining the harmonic SMD force that interacts with the
      specified atom.
    """
    if uses_pbcs:
        smd_force = CustomExternalForce(SMD_FORCE_EXPRESSION_ATOM_PERIODIC)
    else:
        smd_force = CustomExternalForce(SMD_FORCE_EXPRESSION_ATOM_NONPERIODIC)
    smd_force.addGlobalParameter(SMD_FORCE_CONSTANT_PARAMETER_NAME, force_constant)
    smd_force.addPerParticleParameter("x0")
    smd_force.addPerParticleParameter("y0")
    smd_force.addPerParticleParameter("z0")
    smd_force.setForceGroup(31)

    return smd_force


def calculate_com(atom_positions: np.ndarray, atom_masses: np.ndarray) -> np.ndarray:
    """
    Calculate the centre of mass of a group of N atoms, given their positions and masses.

    :param atom_positions: NumPy array of atom positions with dimensions (N, 3)
    :param atom_masses: NumPy array of atomic masses (AMU) with dimension (N)
    :return: NumPy array containing the position of the centre of mass of the atoms with dimension (3)
    """
    # TODO: Make sure this can handle periodic boundary conditions correctly
    #  ^this may not be necessary, as getPositions() should return unwrapped
    #  positions by default!
    return np.sum(
        np.multiply(np.transpose(atom_positions), atom_masses), axis=1
    ) / np.sum(atom_masses)
