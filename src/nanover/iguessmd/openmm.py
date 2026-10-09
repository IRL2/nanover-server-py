import warnings
import os.path
from os import PathLike
from pathlib import Path
from typing import Union, Any
from abc import abstractmethod

import numpy as np

from openmm import CustomExternalForce, CustomCentroidBondForce, OpenMMException
from openmm.app import Simulation

from nanover.openmm import serializer

from nanover.iguessmd.utils import get_every_nth, calculate_unit_tangents, calculate_com

iGUESSMD_FORCE_CONSTANT_PARAMETER_NAME = "smd_k"
iGUESSMD_FORCE_CONSTANT_PARALLEL_PARAMETER_NAME = (
    f"{iGUESSMD_FORCE_CONSTANT_PARAMETER_NAME}_par"
)
iGUESSMD_FORCE_CONSTANT_PERPENDICULAR_PARAMETER_NAME = (
    f"{iGUESSMD_FORCE_CONSTANT_PARAMETER_NAME}_perp"
)
iGUESSMD_FORCE_EXPRESSION_SPHERICAL = (
    f"0.5 * {iGUESSMD_FORCE_CONSTANT_PARAMETER_NAME} * (dx^2 + dy^2 + dz^2)"
)
# tx, ty, tz are components of the unit tangent to the iGUESSMD path at the point x0, y0, z0
iGUESSMD_FORCE_EXPRESSION_PARALLEL = f"0.5 * {iGUESSMD_FORCE_CONSTANT_PARALLEL_PARAMETER_NAME} * (tx*dx + ty*dy + tz*dz)^2"
iGUESSMD_FORCE_EXPRESSION_PERPENDICULAR = f"0.5 * {iGUESSMD_FORCE_CONSTANT_PERPENDICULAR_PARAMETER_NAME} * ((dx^2 + dy^2 + dz^2) - (tx*dx + ty*dy + tz*dz)^2)"
iGUESSMD_FORCE_EXPRESSION_ATOM_NONPERIODIC = f"{iGUESSMD_FORCE_EXPRESSION_PARALLEL} + {iGUESSMD_FORCE_EXPRESSION_PERPENDICULAR}; dx=(x-x0); dy=(y-y0); dz=(z-z0)"
iGUESSMD_FORCE_EXPRESSION_ATOM_PERIODIC = f"{iGUESSMD_FORCE_EXPRESSION_PARALLEL} + {iGUESSMD_FORCE_EXPRESSION_PERPENDICULAR}; dx=(raw_dx - Lx*floor((raw_dx/Lx) + 0.5)); raw_dx=(x-x0); dy=(raw_dy - Ly*floor((raw_dy/Ly) + 0.5)); raw_dy=(y-y0); dz=(raw_dz - Lz*floor((raw_dz/Lz) + 0.5)); raw_dz=(z-z0)"
iGUESSMD_FORCE_EXPRESSION_COM_NONPERIODIC = f"{iGUESSMD_FORCE_EXPRESSION_PARALLEL} + {iGUESSMD_FORCE_EXPRESSION_PERPENDICULAR}; dx=(x1-x0); dy=(y1-y0); dz=(z1-z0)"
iGUESSMD_FORCE_EXPRESSION_COM_PERIODIC = f"{iGUESSMD_FORCE_EXPRESSION_PARALLEL} + {iGUESSMD_FORCE_EXPRESSION_PERPENDICULAR}; dx=(raw_dx - Lx*floor((raw_dx/Lx) + 0.5)); raw_dx=(x1-x0); dy=(raw_dy - Ly*floor((raw_dy/Ly) + 0.5)); raw_dy=(y1-y0); dz=(raw_dz - Lz*floor((raw_dz/Lz) + 0.5)); raw_dz=(z1-z0)"


class OMMiGUESSMDSimulation:
    """
    A wrapper for performing iGUESSMD on an OpenMM simulation.

    This base class defines much of the functionality required for an iGUESSMD
    simulation to be performed using OpenMM, and automatically returns an
    instance of the appropriate subclass, depending on the mode of interaction
    required to perform iGUESSMD (i.e. either single atom or COM).
    """

    @classmethod
    def from_simulation(
        cls,
        simulation: Simulation,
        iguessmd_atom_indices: np.ndarray,
        iguessmd_path: np.ndarray,
        iguessmd_force_constant: float | np.ndarray,
        *,
        name: str | None = None,
    ):
        """
        Construct the iGUESSMD simulation from an existing OpenMM simulation.

        :param simulation: An existing OpenMM Simulation
        :param iguessmd_atom_indices: A NumPy array of indices of the atoms to which the
          iGUESSMD force should be applied (0-D for single atom, 1-D for COM)
        :param iguessmd_path: A NumPy array of coordinates defining the path that the
          iGUESSMD force will take during the iGUESSMD simulation
        :param iguessmd_force_constant: The force constant of the iGUESSMD force
        :param name: An optional name for the simulation instead of default
        """

        # Create instance of iGUESSMD simulation based on type of indices passed
        assert iguessmd_atom_indices.size >= 1
        if iguessmd_atom_indices.size > 1:
            sim = super(cls, OMMiGUESSMDSimulationCOM).__new__(OMMiGUESSMDSimulationCOM)
        else:
            sim = super(cls, OMMiGUESSMDSimulationAtom).__new__(
                OMMiGUESSMDSimulationAtom
            )

        sim.name = name
        sim.simulation = simulation

        # Check if simulation employs periodic boundary conditions
        sim._pbc_box_lengths = None
        sim._sim_uses_pbcs = sim.simulation.system.usesPeriodicBoundaryConditions()
        if sim._sim_uses_pbcs:
            # TODO: only works with orthorhombic PBCs at the moment, need to generalise
            sim._pbc_box_lengths = np.diag(
                sim.simulation.context.getState()
                .getPeriodicBoxVectors(asNumpy=True)
                ._value
            )

        # Initialise all objects relevant to the iGUESSMD simulation
        sim._initialise_iguessmd_simulation(
            iguessmd_atom_indices, iguessmd_path, iguessmd_force_constant
        )

        # Create a checkpoint of the simulation
        sim.checkpoint = sim.simulation.context.createCheckpoint()

        return sim

    @classmethod
    def from_xml_path(
        cls,
        path: PathLike[str],
        iguessmd_atom_indices: np.ndarray,
        iguessmd_path: np.ndarray,
        iguessmd_force_constant: float | np.ndarray,
        *,
        name: str | None = None,
    ):
        """
        Construct the iGUESSMD simulation from an existing NanoVer OpenMM XML file located at a given path.

        :param path: Path of the NanoVer OpenMM XML file
        :param iguessmd_atom_indices: The indices of the atoms to which the iGUESSMD force
          should be applied
        :param iguessmd_path: A NumPy array of coordinates defining the path that the
          iGUESSMD force will take during the iGUESSMD simulation
        :param iguessmd_force_constant: The force constant of the iGUESSMD force
        :param name: An optional name for the simulation instead of default
        """

        if iguessmd_atom_indices.size > 1:
            sim = super(cls, OMMiGUESSMDSimulationCOM).__new__(OMMiGUESSMDSimulationCOM)
        elif iguessmd_atom_indices.size == 1:
            sim = super(cls, OMMiGUESSMDSimulationAtom).__new__(
                OMMiGUESSMDSimulationAtom
            )

        if name is None:
            sim.name = Path(path).stem
        else:
            sim.name = name

        sim.xml_path = path

        # Load the simulation from the path
        with open(sim.xml_path) as infile:
            sim.simulation = serializer.deserialize_simulation(infile)

        # Check if simulation employs periodic boundary conditions
        sim._pbc_box_lengths = None
        sim._sim_uses_pbcs = sim.simulation.system.usesPeriodicBoundaryConditions()
        if sim._sim_uses_pbcs:
            # TODO: only works with orthorhombic PBCs at the moment, need to generalise
            sim._pbc_box_lengths = np.diag(
                sim.simulation.context.getState()
                .getPeriodicBoxVectors(asNumpy=True)
                ._value
            )

        # Initialise all objects relevant to the iGUESSMD simulation
        sim._initialise_iguessmd_simulation(
            iguessmd_atom_indices, iguessmd_path, iguessmd_force_constant
        )

        # Create a checkpoint of the simulation
        sim.checkpoint = sim.simulation.context.createCheckpoint()

        return sim

    def __init__(self, name: str | None = None):
        self.name = name or "Unnamed OpenMM iGUESSMD Simulation"

        self.xml_path: PathLike[str] | None = None

        self.simulation: Simulation | None = None
        self.iguessmd_atom_indices: np.ndarray | None = None
        self.iguessmd_path: np.ndarray | None = None
        self.iguessmd_path_tangents: np.ndarray | None = None
        self.iguessmd_force_constant: float | None = None
        self.iguessmd_force_constant_parallel: float | None = None
        self.iguessmd_force_constant_perpendicular: float | None = None

        self.loaded_iguessmd_force_from_sim: bool = False
        self.n_iguessmd_atom_indices: int | None = None

        self.iguessmd_force: (
            Union[CustomExternalForce, CustomCentroidBondForce] | None
        ) = None

        self.checkpoint: Any | None = None

        self.current_iguessmd_force_position: np.ndarray | None = None
        self.current_iguessmd_force_tangent: np.ndarray | None = None
        self.current_iguessmd_force_position_index: int | None = None
        self.iguessmd_simulation_atom_positions: np.ndarray | None = None
        self.iguessmd_simulation_forces: np.ndarray | None = None
        self.iguessmd_simulation_work_done: np.ndarray | None = None

        self._sim_uses_pbcs: bool | None = None
        self._pbc_box_lengths: np.array | None = None

    @abstractmethod
    def define_iguessmd_simulation_atom_positions_array(self):
        """
        Define the array to which the positions of the atoms with which the
        iGUESSMD force interacts will be saved over the course of the iGUESSMD simulation.
        """
        pass

    @abstractmethod
    def check_for_existing_iguessmd_force(self):
        """
        Check whether the loaded simulation already contains an iGUESSMD force of
        the correct type for the simulation type.
        """
        pass

    @abstractmethod
    def add_iguessmd_force_to_system(self):
        """
        Add the required iGUESSMD force to the system, depending on the type of
        iGUESSMD interaction required (single atom or centre-of-mass).
        """
        pass

    @abstractmethod
    def update_iguessmd_force_position(self):
        """
        Update the position of the iGUESSMD force.
        """
        pass

    @abstractmethod
    def calculate_cumulative_work_done(self):
        """
        Calculate the cumulative work done along the reaction coordinate
        during the iGUESSMD simulation.
        """
        pass

    def _initialise_iguessmd_simulation(
        self,
        iguessmd_atom_indices: np.ndarray,
        iguessmd_path: np.ndarray,
        iguessmd_force_constant: float,
    ):
        """
        Set the fields relevant to the iGUESSMD simulation and add the iGUESSMD force to the system.
        Called upon when constructing an OMMiGUESSMDSimulation using the classmethods
        from_simulation or from_xml_path.
        """
        self.iguessmd_atom_indices = iguessmd_atom_indices
        self.n_iguessmd_atom_indices = self.iguessmd_atom_indices.size
        self.iguessmd_path = iguessmd_path
        if (
            isinstance(iguessmd_force_constant, np.ndarray)
            and iguessmd_force_constant.size == 2
        ):
            self.iguessmd_force_constant = iguessmd_force_constant
            self.iguessmd_force_constant_parallel = iguessmd_force_constant[0]
            self.iguessmd_force_constant_perpendicular = iguessmd_force_constant[1]
        else:
            self.iguessmd_force_constant = self.iguessmd_force_constant_parallel = (
                self.iguessmd_force_constant_perpendicular
            ) = iguessmd_force_constant

        self.define_iguessmd_simulation_atom_positions_array()

        self.current_iguessmd_force_position = self.iguessmd_path[0]
        self.current_iguessmd_force_position_index = 0

        self.calculate_iguessmd_path_tangents()
        self.current_iguessmd_force_tangent = self.iguessmd_path_tangents[0]

        # Check whether iGUESSMD force is already present
        self.check_for_existing_iguessmd_force()

        if not self.loaded_iguessmd_force_from_sim:
            # Create iGUESSMD force and add it to the system
            self.add_iguessmd_force_to_system()

        self.simulation.context.reinitialize(preserveState=True)

    def reset(self):
        """
        Reset the iGUESSMD simulation to its initial state, and reset the arrays output
        by the iGUESSMD simulation.
        """
        assert (
            self.simulation is not None
            and self.checkpoint is not None
            and self.iguessmd_path is not None
            and self.iguessmd_force is not None
        )

        # Reset iGUESSMD force position
        self.current_iguessmd_force_position = self.iguessmd_path[0]
        self.current_iguessmd_force_position_index = 0
        self.update_iguessmd_force_position()

        # Reset simulation, reinitialise context to be safe
        self.simulation.context.reinitialize()
        self.simulation.context.loadCheckpoint(self.checkpoint)

        # Reset or remove iGUESSMD simulation arrays
        self.define_iguessmd_simulation_atom_positions_array()
        del self.iguessmd_simulation_forces
        del self.iguessmd_simulation_work_done

    def remove_iguessmd_force_from_system(self):
        """
        Remove any iGUESSMD forces from the system.
        """
        forces = self.simulation.system.getForces()
        forces_to_remove = []
        for i in range(len(forces)):
            if type(forces[i]) is type(self.iguessmd_force):
                try:
                    if (
                        forces[i].getGlobalParameterName(0)
                        == iGUESSMD_FORCE_CONSTANT_PARALLEL_PARAMETER_NAME
                    ):
                        forces_to_remove.append(i)
                except OpenMMException:
                    continue
                # forces_to_remove.append(i)

        # Remove any iGUESSMD forces, accounting for the changes in indices
        # as forces are removed
        for j in range(len(forces_to_remove)):
            self.simulation.system.removeForce(forces_to_remove[j] - j)

    def calculate_iguessmd_path_tangents(self):
        """
        Calculate the unit tangent vectors of the iGUESSMD path using the forward
        difference approximation.
        """

        # Check that path exists
        assert self.iguessmd_path is not None

        # Calculate unit tangents
        self.iguessmd_path_tangents = calculate_unit_tangents(self.iguessmd_path)

    def get_iguessmd_atom_positions(self):
        """
        Retrieve the positions of the atoms with which the iGUESSMD force is
        interacting, and add them to the array of positions to save.
        """
        positions = self.simulation.context.getState(
            getPositions=True, enforcePeriodicBox=False
        ).getPositions(asNumpy=True)
        self.iguessmd_simulation_atom_positions[
            self.current_iguessmd_force_position_index
        ] = positions[self.iguessmd_atom_indices]

    def run_equilibration_with_initial_restraint(self, n_steps: int):
        """
        Perform an equilibration of the system with the restraint fixed at the
        initial position defined in the iGUESSMD force path.

        :param n_steps: Number of simulation steps to run equilibration for.
        """

        try:
            assert np.all(self.current_iguessmd_force_position == self.iguessmd_path[0])
        except AssertionError:
            raise AssertionError(
                "Restraint is not located at the initial position "
                "of the iGUESSMD force, equilibration aborted. To prepare the system "
                "appropriately, the restraint should be placed at the first point "
                "on the path that the iGUESSMD force will take during the iGUESSMD simulation."
            )

        self.simulation.step(n_steps)

    def save_simulation(
        self,
        output_filepath: PathLike[str],
        save_state: bool = False,
        save_iguessmd_force: bool | None = False,
    ):
        """
        Save the simulation to a NanoVer OpenMM XML file, with the option to include the
        iGUESSMD force in the XML file.

        :param output_filepath: Path to output file to save the simulation to.
        :param save_state: If True, save the present state of the simulation to the XML file.
        :param save_iguessmd_force: Bool defining whether to save the iGUESSMD force in the XML file (Optional).
        """
        assert output_filepath is not None

        if save_iguessmd_force:
            with open(output_filepath, "w") as outfile:
                outfile.write(
                    serializer.serialize_simulation(
                        self.simulation, save_state=save_state
                    )
                )

        else:
            # Temporarily remove iGUESSMD force from simulation
            self.remove_iguessmd_force_from_system()

            # Save simulation without SMD force
            with open(output_filepath, "w") as outfile:
                xml_string = serializer.serialize_simulation(
                    self.simulation, save_state=save_state
                )
                # Manually remove parameters for now...
                xml_string_lines = xml_string.splitlines()
                for line in xml_string_lines:
                    # Remove only parameter left over from iGUESSMD force
                    if (
                        iGUESSMD_FORCE_CONSTANT_PARAMETER_NAME in line
                        or iGUESSMD_FORCE_CONSTANT_PARALLEL_PARAMETER_NAME in line
                        or iGUESSMD_FORCE_CONSTANT_PERPENDICULAR_PARAMETER_NAME in line
                    ):
                        line_index = xml_string_lines.index(line)
                        n_tabs = line.index("<")
                        param_string = line.split("/")
                        param_string_list = sum(
                            [string.split() for string in param_string], []
                        )
                        for substring in param_string_list:
                            if iGUESSMD_FORCE_CONSTANT_PARAMETER_NAME in substring:
                                param_string_list.remove(substring)
                        # TODO: make more efficient
                        for substring in param_string_list:
                            if (
                                iGUESSMD_FORCE_CONSTANT_PARALLEL_PARAMETER_NAME
                                in substring
                            ):
                                param_string_list.remove(substring)
                        for substring in param_string_list:
                            if (
                                iGUESSMD_FORCE_CONSTANT_PERPENDICULAR_PARAMETER_NAME
                                in substring
                            ):
                                param_string_list.remove(substring)
                        if self._sim_uses_pbcs:
                            for substring in param_string_list:
                                if "Lx" in substring:
                                    param_string_list.remove(substring)
                            for substring in param_string_list:
                                if "Ly" in substring:
                                    param_string_list.remove(substring)
                            for substring in param_string_list:
                                if "Lz" in substring:
                                    param_string_list.remove(substring)
                        xml_string_lines[line_index] = (
                            n_tabs * "\t"
                            + " ".join(param_string_list[:-1])
                            + "/"
                            + param_string_list[-1]
                        )
                xml_string = "\n".join(xml_string_lines)
                outfile.write(xml_string)

            # Add iGUESSMD force back to the system
            self.add_iguessmd_force_to_system()

    def generate_starting_structures(
        self,
        interval_ps: float,
        n_structures: int,
        output_directory: PathLike[str] | None = None,
        filename_prefix: str | None = None,
        save_iguessmd_force: bool | None = None,
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
        :param save_iguessmd_force: Bool defining whether to save the iGUESSMD force in the XML file (Optional).
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
                save_iguessmd_force=save_iguessmd_force,
            )

        print(f"Structure generation complete: {n_structures} structures generated.")

    def run_iguessmd(self, progress_interval: int | None = 100000):
        """
        Perform an iGUESSMD simulation on the system using the iGUESSMD force and
        path defined, store the positions of the atoms with which the
        iGUESSMD force interacts, and calculate the work done along the reaction
        coordinate.

        :param progress_interval: Interval defining how regularly to print the progress
          of the simulation, in terms of the number of simulation steps.
        """

        # Get initial atom positions for iGUESSMD force at initial position
        assert (
            self.iguessmd_simulation_atom_positions is not None
            and self.iguessmd_force is not None
            and np.all(self.current_iguessmd_force_position == self.iguessmd_path[0])
            and self.current_iguessmd_force_position_index == 0
        )
        self.get_iguessmd_atom_positions()

        # Run iGUESSMD procedure
        n_steps = self.iguessmd_path.shape[0]

        print("Starting SMD simulation...")
        for step in range(1, n_steps):
            # Perform single simulation step
            self.simulation.step(1)

            # Update force position index
            self.current_iguessmd_force_position_index = step

            # Update iGUESSMD force position
            self.update_iguessmd_force_position()

            # Retrieve atom positions on step
            self.get_iguessmd_atom_positions()

            # Print step at intervals
            if step % progress_interval == 0:
                print(f"Steps completed: {step}")

        print("iGUESSMD simulation completed. Calculating work done...")

        # Calculate the work done along the reaction coordinate
        self.calculate_cumulative_work_done()

        print("Work done calculated.")

    def _calculate_iguessmd_forces(self, interaction_centre_positions):
        """
        Calculate the iGUESSMD forces that acted on the system during the simulation in kJ mol-1 nm-1.

        :param interaction_centre_positions: Array of positions defining the centre
          (single atom or COM of group of atoms) with which the iGUESSMD force interacted during the simulation.
        """

        assert np.all(self.iguessmd_path.shape == interaction_centre_positions.shape)

        displacements = interaction_centre_positions - self.iguessmd_path

        # Calculate force component along RC
        parallel_forces = (
            -self.iguessmd_force_constant_parallel
            * (
                np.linalg.vecdot(displacements, self.iguessmd_path_tangents, axis=1)
                / np.linalg.vecdot(
                    self.iguessmd_path_tangents, self.iguessmd_path_tangents, axis=1
                )
            ).reshape((self.iguessmd_path_tangents.shape[0], 1))
            * self.iguessmd_path_tangents
        )
        perpendicular_forces = -self.iguessmd_force_constant_perpendicular * (
            displacements
            - (
                (
                    np.linalg.vecdot(displacements, self.iguessmd_path_tangents, axis=1)
                    / np.linalg.vecdot(
                        self.iguessmd_path_tangents, self.iguessmd_path_tangents, axis=1
                    )
                ).reshape((self.iguessmd_path_tangents.shape[0], 1))
                * self.iguessmd_path_tangents
            )
        )

        self.iguessmd_simulation_forces = parallel_forces + perpendicular_forces

    def _calculate_work_done(self):
        """
        Calculate the work done along the reaction coordinate in kJ mol-1.
        """
        assert (
            self.iguessmd_path is not None
            and self.iguessmd_simulation_forces is not None
        )
        smd_force_displacements = np.diff(self.iguessmd_path, axis=0)
        work_done_array = np.zeros(self.iguessmd_simulation_forces.shape[0])
        for i in range(smd_force_displacements.shape[0]):
            work_done_array[i + 1] = np.dot(
                self.iguessmd_simulation_forces[i], smd_force_displacements[i]
            )
        self.iguessmd_simulation_work_done = np.cumsum(work_done_array, axis=0)

    def save_iguessmd_simulation_data(
        self,
        path: PathLike[str] = None,
        every_nth: int | None = None,
        include_end: bool = False,
        save_work_done: bool = True,
        save_atom_positions: bool = True,
        work_done_dtype: np.dtype = np.float32,
        atom_positions_dtype: np.dtype = np.float32,
    ):
        """
        Save the data produced by the iGUESSMD simulation in binary form that can be read
        into NumPy arrays. The following data are saved, in the order listed below:

        - Trajectories of the atoms to which the iGUESSMD force was applied, in nm
        - Work done along the reaction coordinate defined by the path of the iGUESSMD force, in kJ mol-1

        :param path: Path to the file to which the data will be saved.
        :param every_nth:  (int | None) only save the values of work and/or positions at every
          nth point along the trajectory.
        :param include_end: (Bool) whether to save final values of the work and/or positions
          regardless of stride defined by every nth
        :param save_work_done: Bool determining whether to save the work done
        :param save_atom_positions: Bool determining whether to save the positions of the atom(s)
        :param work_done_dtype: Data type of the work done array to save.
        :param atom_positions_dtype: Data type of the atom positions array to save.
        """

        if path is None:
            raise ValueError(
                "Output file path cannot be None. Please specify an output file path."
            )

        elif self.iguessmd_simulation_work_done is None:
            raise ValueError(
                "Missing values for the work done. This data can only be saved after "
                "the iGUESSMD calculation is completed."
            )

        elif self.iguessmd_simulation_atom_positions is None or np.all(
            self.iguessmd_simulation_atom_positions == 0.0
        ):
            raise ValueError(
                "Missing values for the atom positions. This data can only be saved after "
                "the iGUESSMD calculation is completed."
            )

        if every_nth is None:
            every_nth = 1

        if (
            every_nth is not None
            and include_end
            and (self.iguessmd_simulation_work_done.size - 1) % every_nth != 0
        ):
            warnings.warn(
                "Choice of every_nth yields different time step between the "
                "final two array entries compared to the rest of the trajectory."
            )


        iguessmd_simulation_data = {}
        iguessmd_simulation_data["data_timestep_ps"] = every_nth * self.simulation.integrator.getStepSize()._value

        # Optionally save atom positions
        if save_atom_positions:
            iguessmd_simulation_data["iguessmd_simulation_atom_positions"] = get_every_nth(
                self.iguessmd_simulation_atom_positions,
                axis=0,
                every_nth=every_nth,
                include_end=include_end,
            ).astype(atom_positions_dtype)
            print("Atom positions saved to simulation data file.")

        # Optionally save work done
        if save_work_done:
            iguessmd_simulation_data["iguessmd_simulation_work_done"] = get_every_nth(
                self.iguessmd_simulation_work_done.astype(work_done_dtype),
                axis=0,
                every_nth=every_nth,
                include_end=include_end,
            )
            print("Work done saved to simulation data file.")

        np.savez_compressed(path, **iguessmd_simulation_data)


    def save_general_iguessmd_data(self, path: PathLike | str = None):
        """
        Saves general data related to the iGUESSMD simulation in binary form that can be read
        into NumPy arrays. The following data are saved, in the order listed below:

        - Indices of the atoms to which the iGUESSMD force is applied
        - Positions defining the path of the iGUESSMD force, in nm
        - Force constant of the iGUESSMD force, in kJ mol-1 nm-2
        - Temperature of the simulation, in Kelvin
        - Time step of the simulation, in picoseconds

        :param path: Path to the file to which the data will be saved.
        """

        if path is None:
            raise ValueError(
                "Output file path cannot be None. Please specify an output file path."
            )

        # Create dictionary containing general iGUESSMD data
        general_iguessmd_data = {
            "iguessmd_atom_indices" : self.iguessmd_atom_indices,
            "iguessmd_path" : self.iguessmd_path,
            "iguessmd_force_constant" : self.iguessmd_force_constant,
            "temperature_K" : self.simulation.integrator.getTemperature()._value,
            "timestep_ps" : self.simulation.integrator.getStepSize()._value,
        }

        # Save dictionary to file defined by path
        np.savez_compressed(path, **general_iguessmd_data)


class OMMiGUESSMDSimulationAtom(OMMiGUESSMDSimulation):
    """
    A class for performing constant velocity iGUESSMD on an OpenMM simulations,
    where the iGUESSMD force is applied to a single atom.
    """

    def __init__(self, name: str | None = None):
        super().__init__(name)

    def check_for_existing_iguessmd_force(self):
        try:
            force_constant_parallel = self.simulation.context.getParameter(
                iGUESSMD_FORCE_CONSTANT_PARALLEL_PARAMETER_NAME
            )
            force_constant_perpendicular = self.simulation.context.getParameter(
                iGUESSMD_FORCE_CONSTANT_PERPENDICULAR_PARAMETER_NAME
            )
            n_forces = self.simulation.system.getNumForces()
            iguessmd_force = self.simulation.system.getForce(n_forces - 1)
            params = iguessmd_force.getParticleParameters(0)
            assert isinstance(iguessmd_force, CustomExternalForce)
            assert (
                iguessmd_force.getEnergyFunction()
                == iGUESSMD_FORCE_EXPRESSION_ATOM_PERIODIC
                or iguessmd_force.getEnergyFunction()
                == iGUESSMD_FORCE_EXPRESSION_ATOM_NONPERIODIC
            )
            assert iguessmd_force.getNumParticles() == 1
            assert force_constant_parallel == self.iguessmd_force_constant_parallel
            assert (
                force_constant_perpendicular
                == self.iguessmd_force_constant_perpendicular
            )
            assert params[0] == self.iguessmd_atom_indices
            assert np.all(
                params[1]
                == np.array(
                    [
                        *self.current_iguessmd_force_position,
                        *self.current_iguessmd_force_tangent,
                    ]
                )
            )
            print("iGUESSMD force already present in loaded simulation.")
            self.iguessmd_force = iguessmd_force
            self.loaded_iguessmd_force_from_sim = True

        except OpenMMException:
            self.loaded_iguessmd_force_from_sim = False

    def add_iguessmd_force_to_system(self):
        """
        Add an iGUESSMD force to the OpenMM system that interacts with the
        specified atom.
        """

        x0, y0, z0 = self.current_iguessmd_force_position
        tx, ty, tz = self.current_iguessmd_force_tangent
        iguessmd_force = iguessmd_single_atom_force(
            self.iguessmd_force_constant_parallel,
            self.iguessmd_force_constant_perpendicular,
            self._sim_uses_pbcs,
            self._pbc_box_lengths,
        )
        iguessmd_force.addParticle(self.iguessmd_atom_indices, [x0, y0, z0, tx, ty, tz])
        self.iguessmd_force = iguessmd_force
        self.simulation.system.addForce(self.iguessmd_force)
        self.simulation.context.reinitialize(preserveState=True)

    def update_iguessmd_force_position(self):
        self.current_iguessmd_force_position = self.iguessmd_path[
            self.current_iguessmd_force_position_index
        ]
        self.current_iguessmd_force_tangent = self.iguessmd_path_tangents[
            self.current_iguessmd_force_position_index
        ]

        x0, y0, z0 = self.current_iguessmd_force_position
        tx, ty, tz = self.current_iguessmd_force_tangent
        self.iguessmd_force.setParticleParameters(
            0, self.iguessmd_atom_indices, [x0, y0, z0, tx, ty, tz]
        )
        self.iguessmd_force.updateParametersInContext(self.simulation.context)

    def define_iguessmd_simulation_atom_positions_array(self):
        self.iguessmd_simulation_atom_positions = np.zeros(
            (self.iguessmd_path.shape[0], 3)
        )

    def calculate_cumulative_work_done(self):
        """
        Calculate the cumulative work done by the iGUESSMD force on the
        atom with which it interacts over the iGUESSMD simulation.
        """
        assert not np.array_equal(
            self.iguessmd_simulation_atom_positions,
            np.zeros((self.iguessmd_path.shape[0], 3)),
        )
        self._calculate_iguessmd_forces(self.iguessmd_simulation_atom_positions)
        self._calculate_work_done()


class OMMiGUESSMDSimulationCOM(OMMiGUESSMDSimulation):
    """
    A class for performing constant velocity iGUESSMD on an OpenMM simulations,
    where the iGUESSMD force is applied to the centre of mass of a specified
    group of atoms.
    """

    def __init__(self, name: str | None = None):
        super().__init__(name)
        self.iguessmd_com_positions: np.ndarray | None = None

    def check_for_existing_iguessmd_force(self):
        try:
            force_constant_parallel = self.simulation.context.getParameter(
                iGUESSMD_FORCE_CONSTANT_PARALLEL_PARAMETER_NAME
            )
            force_constant_perpendicular = self.simulation.context.getParameter(
                iGUESSMD_FORCE_CONSTANT_PERPENDICULAR_PARAMETER_NAME
            )
            n_forces = self.simulation.system.getNumForces()
            iguessmd_force = self.simulation.system.getForce(n_forces - 1)
            assert isinstance(iguessmd_force, CustomCentroidBondForce)
            assert (
                iguessmd_force.getEnergyFunction()
                == iGUESSMD_FORCE_EXPRESSION_COM_NONPERIODIC
                or iguessmd_force.getEnergyFunction()
                == iGUESSMD_FORCE_EXPRESSION_COM_PERIODIC
            )
            assert iguessmd_force.getNumGroups() == 1
            assert force_constant_parallel == self.iguessmd_force_constant_parallel
            assert (
                force_constant_perpendicular
                == self.iguessmd_force_constant_perpendicular
            )
            assert np.all(
                iguessmd_force.getGroupParameters(0)[0] == self.iguessmd_atom_indices
            )
            assert np.all(
                iguessmd_force.getBondParameters(0)[1]
                == np.array(
                    [
                        *self.current_iguessmd_force_position,
                        *self.current_iguessmd_force_tangent,
                    ]
                )
            )
            print("iGUESSMD force already present in loaded simulation.")
            self.iguessmd_force = iguessmd_force
            self.loaded_iguessmd_force_from_sim = True

        except OpenMMException:
            self.loaded_iguessmd_force_from_sim = False

    def add_iguessmd_force_to_system(self):
        """
        Add an iGUESSMD force to the OpenMM system that interacts with the
        centre of mass of the specified group of atoms.
        """

        x0, y0, z0 = self.current_iguessmd_force_position
        tx, ty, tz = self.current_iguessmd_force_tangent
        iguessmd_force = iguessmd_com_force(
            self.iguessmd_force_constant_parallel,
            self.iguessmd_force_constant_perpendicular,
            self._sim_uses_pbcs,
            self._pbc_box_lengths,
        )
        iguessmd_force.addGroup(self.iguessmd_atom_indices)
        iguessmd_force.addBond([0], [x0, y0, z0, tx, ty, tz])
        self.iguessmd_force = iguessmd_force
        self.simulation.system.addForce(self.iguessmd_force)
        self.simulation.context.reinitialize(preserveState=True)

    def update_iguessmd_force_position(self):
        self.current_iguessmd_force_position = self.iguessmd_path[
            self.current_iguessmd_force_position_index
        ]
        self.current_iguessmd_force_tangent = self.iguessmd_path_tangents[
            self.current_iguessmd_force_position_index
        ]
        x0, y0, z0 = self.current_iguessmd_force_position
        tx, ty, tz = self.current_iguessmd_force_tangent
        self.iguessmd_force.setBondParameters(0, [0], [x0, y0, z0, tx, ty, tz])
        self.iguessmd_force.updateParametersInContext(self.simulation.context)

    def define_iguessmd_simulation_atom_positions_array(self):
        self.iguessmd_simulation_atom_positions = np.zeros(
            (self.iguessmd_path.shape[0], self.iguessmd_atom_indices.size, 3)
        )

    def _calculate_com(self, atom_positions: np.ndarray, atom_masses: np.ndarray):
        """
        Calculate the centre of mass of the N atoms to which the iGUESSMD force has been applied.

        :param atom_positions: NumPy array of atom positions with dimensions (N, 3)
        :param atom_masses: NumPy array of atomic masses (AMU) with dimension (N)
        :return: NumPy array containing the position of the centre of mass of the atoms with dimension (3)
        """
        assert np.all(
            atom_positions.shape == np.array((self.iguessmd_atom_indices.size, 3))
        )
        return calculate_com(atom_positions, atom_masses)

    def _calculate_com_trajectory(self):
        """
        Calculate the trajectory that the COM follows during the iGUESSMD simulation.
        """
        assert not np.array_equal(
            self.iguessmd_simulation_atom_positions,
            np.zeros((self.iguessmd_path.shape[0], self.iguessmd_atom_indices.size, 3)),
        )
        atom_masses = np.zeros(self.n_iguessmd_atom_indices)
        for index in range(self.n_iguessmd_atom_indices):
            atom_masses[index] = self.simulation.system.getParticleMass(
                self.iguessmd_atom_indices[index]
            )._value
        self.iguessmd_com_positions = np.array(
            [
                self._calculate_com(
                    self.iguessmd_simulation_atom_positions[i], atom_masses
                )
                for i in range(self.iguessmd_path.shape[0])
            ]
        )

    def calculate_cumulative_work_done(self):
        """
        Calculate the cumulative work done by the iGUESSMD force on the
        COM of the atoms with which it interacts over the iGUESSMD simulation.
        """
        self._calculate_com_trajectory()
        self._calculate_iguessmd_forces(self.iguessmd_com_positions)
        self._calculate_work_done()

    def save_iguessmd_simulation_data(
        self,
        path: PathLike[str] = None,
        every_nth: int | None = None,
        include_end: bool = False,
        save_com_positions: bool = True,
        com_positions_dtype: np.dtype = np.float32,
        **kwargs,
    ):
        """
        Save the data produced by the iGUESSMD simulation in binary form that can be read
        into NumPy arrays. The following data can be saved (optionally), in the order listed below:

        - Trajectories of the atoms to which the iGUESSMD force was applied, in nm
        - Work done along the reaction coordinate defined by the path of the iGUESSMD force, in kJ mol-1
        - Trajectories of the COM of the atoms to which the iGUESSMD force was applied, in nm

        :param path: Path to the file to which the data will be saved.
        :param every_nth:  (int | None) only save the values of work and/or positions at every
          nth point along the trajectory.
        :param include_end: (Bool) whether to save final values of the work and/or positions
          regardless of stride defined by every nth
        :param save_work_done: Bool determining whether to save the work done
        :param save_atom_positions: Bool determining whether to save the positions of the atom(s)
        :param save_com_positions: Bool determining whether to save the positions of the COM
        :param work_done_dtype: Data type of the work done array to save.
        :param atom_positions_dtype: Data type of the atom positions array to save.
        :param com_positions_dtype: Data type of the COM positions array to save.
        """
        # TODO: Think about whether there is a cleaner way to achieve this

        # Optionally save work done and atomic coordinates
        super().save_iguessmd_simulation_data(path, every_nth=every_nth, include_end=include_end, **kwargs)

        if self.iguessmd_com_positions is None or np.all(self.iguessmd_com_positions == 0.0):
            raise ValueError(
                "Missing values for the atom positions. This data can only be saved after"
                "the iGUESSMD calculation is completed."
            )

        if every_nth is None:
            every_nth = 1

        if save_com_positions:
            # Load existing saved data
            iguessmd_simulation_data = dict(np.load(path))
            iguessmd_simulation_data["iguessmd_com_positions"] = get_every_nth(
                self.iguessmd_com_positions.astype(com_positions_dtype),
                axis=0,
                every_nth=every_nth,
                include_end=include_end,
            )
            np.savez_compressed(path, **iguessmd_simulation_data)

        # Save

        # # Optionally save COM positions
        # with open(path, "ab+") as outfile:
        #     if save_com_positions:
        #         np.save(
        #             outfile,
        #             get_every_nth(
        #                 self.iguessmd_com_positions.astype(com_positions_dtype),
        #                 axis=0,
        #                 every_nth=every_nth,
        #                 include_end=include_end,
        #             ),
        #         )
        #         print("COM positions saved to simulation data file.")


def iguessmd_com_force(
    parallel_force_constant: float,
    perpendicular_force_constant: float,
    uses_pbcs: bool,
    pbc_box_lengths: np.ndarray | None = None,
):
    """
    Defines a harmonic restraint force for the COM of a group of atoms for performing iGUESSMD.

    :param parallel_force_constant: Force constant of the harmonic restraint to be applied to the
      COM of the group of atoms in the direction tangential to the RC, in units kJ mol-1 nm-2
    :param perpendicular_force_constant: Force constant of the harmonic restraint to be applied to the
      COM of the group of atoms in the direction perpendicular to the RC, in units kJ mol-1 nm-2
    :param uses_pbcs: Bool specifying whether to use periodic boundary conditions for the
      harmonic restraint
    :return: CustomCentroidBondForce defining the harmonic iGUESSMD force that interacts with
      the COM of the specified atoms.
    """

    if uses_pbcs:
        iguessmd_force = CustomCentroidBondForce(
            1, iGUESSMD_FORCE_EXPRESSION_COM_PERIODIC
        )
    else:
        iguessmd_force = CustomCentroidBondForce(
            1, iGUESSMD_FORCE_EXPRESSION_COM_NONPERIODIC
        )

    iguessmd_force.addGlobalParameter(
        f"{iGUESSMD_FORCE_CONSTANT_PARALLEL_PARAMETER_NAME}", parallel_force_constant
    )
    iguessmd_force.addGlobalParameter(
        f"{iGUESSMD_FORCE_CONSTANT_PERPENDICULAR_PARAMETER_NAME}",
        perpendicular_force_constant,
    )
    if uses_pbcs and pbc_box_lengths is not None:
        iguessmd_force.addGlobalParameter("Lx", pbc_box_lengths[0])
        iguessmd_force.addGlobalParameter("Ly", pbc_box_lengths[1])
        iguessmd_force.addGlobalParameter("Lz", pbc_box_lengths[2])
    iguessmd_force.addPerBondParameter("x0")
    iguessmd_force.addPerBondParameter("y0")
    iguessmd_force.addPerBondParameter("z0")
    iguessmd_force.addPerBondParameter("tx")
    iguessmd_force.addPerBondParameter("ty")
    iguessmd_force.addPerBondParameter("tz")
    iguessmd_force.setUsesPeriodicBoundaryConditions(uses_pbcs)
    iguessmd_force.setForceGroup(31)

    return iguessmd_force


def iguessmd_single_atom_force(
    parallel_force_constant: float,
    perpendicular_force_constant: float,
    uses_pbcs: bool,
    pbc_box_lengths: np.ndarray | None = None,
):
    """
    Defines a harmonic restraint force for a single atom for performing iGUESSMD.

    :param parallel_force_constant: Force constant of the harmonic restraint to be applied to the
      specified atom in the direction tangential to the RC, in units kJ mol-1 nm-2
    :param perpendicular_force_constant: Force constant of the harmonic restraint to be applied to the
      specified atom in the direction perpendicular to the RC, in units kJ mol-1 nm-2
    :param uses_pbcs: Bool specifying whether to use periodic boundary conditions for the
      harmonic restraint
    :return: CustomExternalForce defining the harmonic iGUESSMD force that interacts with the
      specified atom.
    """
    if uses_pbcs:
        iguessmd_force = CustomExternalForce(iGUESSMD_FORCE_EXPRESSION_ATOM_PERIODIC)
    else:
        iguessmd_force = CustomExternalForce(iGUESSMD_FORCE_EXPRESSION_ATOM_NONPERIODIC)
    iguessmd_force.addGlobalParameter(
        f"{iGUESSMD_FORCE_CONSTANT_PARALLEL_PARAMETER_NAME}", parallel_force_constant
    )
    iguessmd_force.addGlobalParameter(
        f"{iGUESSMD_FORCE_CONSTANT_PERPENDICULAR_PARAMETER_NAME}",
        perpendicular_force_constant,
    )
    if uses_pbcs and pbc_box_lengths is not None:
        iguessmd_force.addGlobalParameter("Lx", pbc_box_lengths[0])
        iguessmd_force.addGlobalParameter("Ly", pbc_box_lengths[1])
        iguessmd_force.addGlobalParameter("Lz", pbc_box_lengths[2])
    iguessmd_force.addPerParticleParameter("x0")
    iguessmd_force.addPerParticleParameter("y0")
    iguessmd_force.addPerParticleParameter("z0")
    iguessmd_force.addPerParticleParameter("tx")
    iguessmd_force.addPerParticleParameter("ty")
    iguessmd_force.addPerParticleParameter("tz")
    iguessmd_force.setForceGroup(31)

    return iguessmd_force
