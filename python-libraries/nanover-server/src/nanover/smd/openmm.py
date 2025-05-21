import os.path
from os import PathLike
from pathlib import Path
from typing import Optional, Union, Any
from abc import abstractmethod

import numpy as np

from openmm import CustomExternalForce, CustomCentroidBondForce, OpenMMException
from openmm.app import Simulation

from nanover.openmm import serializer


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
        name: Optional[str] = None,
    ):
        """
        Construct the SMD simulation from an existing OpenMM simulation.

        :param simulation: An existing OpenMM Simulation
        :param smd_atom_indices: The indices of the atoms to which the SMD force
          should be applied
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

        sim.smd_atom_indices = smd_atom_indices
        sim.smd_path = smd_path
        sim.smd_force_constant = smd_force_constant
        sim.n_smd_atom_indices = sim.smd_atom_indices.size
        sim.define_smd_simulation_atom_positions_array()

        # Create SMD force and add it to the system
        sim.current_smd_force_position = sim.smd_path[0]
        sim.current_smd_force_position_index = 0

        # Check whether SMD force is already present
        sim.check_for_existing_smd_force()

        if not sim.loaded_smd_force_from_sim:
            # Create SMD force and add it to the system
            sim.add_smd_force_to_system()

        sim.simulation.context.reinitialize(preserveState=True)

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
        name: Optional[str] = None,
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

        sim.smd_atom_indices = smd_atom_indices
        sim.n_smd_atom_indices = sim.smd_atom_indices.size
        sim.smd_path = smd_path
        sim.smd_force_constant = smd_force_constant
        sim.define_smd_simulation_atom_positions_array()

        sim.current_smd_force_position = sim.smd_path[0]
        sim.current_smd_force_position_index = 0

        # Check whether SMD force is already present
        sim.check_for_existing_smd_force()

        if not sim.loaded_smd_force_from_sim:
            # Create SMD force and add it to the system
            sim.add_smd_force_to_system()

        sim.simulation.context.reinitialize(preserveState=True)

        # Create a checkpoint of the simulation
        sim.checkpoint = sim.simulation.context.createCheckpoint()

        return sim

    def __init__(self, name: Optional[str] = None):
        self.name = name or "Unnamed OpenMM SMD Simulation"

        self.xml_path: Optional[PathLike[str]] = None

        self.simulation: Optional[Simulation] = None
        self.smd_atom_indices: Optional[np.ndarray] = None
        self.smd_path: Optional[np.ndarray] = None
        self.smd_force_constant: Optional[float] = None

        self.loaded_smd_force_from_sim: bool = False
        self.n_smd_atom_indices: Optional[int] = None

        self.smd_force: Optional[
            Union[CustomExternalForce, CustomCentroidBondForce]
        ] = None

        self.checkpoint: Optional[Any] = None

        self.current_smd_force_position: Optional[np.ndarray] = None
        self.current_smd_force_position_index: Optional[int] = None
        self.smd_simulation_atom_positions: Optional[np.ndarray] = None
        self.smd_simulation_forces: Optional[np.ndarray] = None
        self.smd_simulation_work_done: Optional[np.ndarray] = None

        self._sim_uses_pbcs: Optional[bool] = None

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

    def reset(self):
        """
        Reset the SMD simulation to its initial state, and reset the arrays output
        by the SMD simulation.
        """
        assert self.simulation is not None and self.checkpoint is not None and self.smd_path is not None and self.smd_force is not None

        # Reset SMD force position
        self.current_smd_force_position = self.smd_path[0]
        self.current_smd_force_position_index = 0
        self.update_smd_force_position()

        # Reset simulation, reinitialise context to be safe
        self.simulation.context.reinitialize()
        self.simulation.context.loadCheckpoint(self.checkpoint)

        # Reset SMD simulation arrays
        self.smd_simulation_atom_positions = None
        self.define_smd_simulation_atom_positions_array()
        self.smd_simulation_forces = None
        self.smd_simulation_work_done = None



    def remove_smd_force_from_system(self):
        """
        Remove any SMD forces from the system.
        """
        forces = self.simulation.system.getForces()
        forces_to_remove = []
        for i in range(len(forces)):
            if (
                type(forces[i]) == type(self.smd_force)
                and forces[i].getGlobalParameterName(0) == "smd_k"
            ):
                forces_to_remove.append(i)

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
        save_smd_force: Optional[bool] = False,
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
                # TODO: Work out how to deal with these unwanted parameters...
                #  WARNING: this will not work if you have other parameters!
                #  maybe change smd_k to a per bond parameter?
                xml_string = "\n".join(
                    x for x in xml_string.splitlines() if "smd_k" not in x
                )
                outfile.write(xml_string)

            # Add SMD force back to the system
            self.add_smd_force_to_system()

    def generate_starting_structures(
        self,
        interval_ps: float,
        n_structures: int,
        output_directory: Optional[PathLike[str]] = None,
        filename_prefix: Optional[str] = None,
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
                output_filepath=outfile_path, save_state=True, save_smd_force=False
            )

        print(f"Structure generation complete: {n_structures} structures generated.")

    def run_smd(self, progress_interval: Optional[int] = 100000):
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
            and np.all(self.current_smd_force_position == self.smd_path[0])
            and self.current_smd_force_position_index == 0
        )
        self.get_smd_atom_positions()

        # Run SMD procedure
        n_steps = self.smd_path.shape[0]

        print("Starting SMD simulation...")
        for step in range(1, n_steps):

            # Update force position index
            self.current_smd_force_position_index = step

            # Update SMD force position
            self.update_smd_force_position()

            # Perform single simulation step
            self.simulation.step(1)

            # Retrieve atom positions on step
            self.get_smd_atom_positions()

            # Print every 10000 steps
            if step % progress_interval == 0:
                print(f"Steps completed: {step}")

        print("SMD simulation completed. Calculating work done...")

        # Calculate the work done along the reaction coordinate
        self.calculate_cumulative_work_done()

        print("Work done calculated.")

    def _calculate_forces(self, interaction_centre_positions):
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

    def save_smd_simulation_data(self, path: PathLike[str] = None):
        """
        Save the data produced by the SMD simulation in binary form that can be read
        into NumPy arrays. The following data are saved, in the order listed below:

        - Trajectories of the atoms to which the SMD force was applied, in nm
        - Work done along the reaction coordinate defined by the path of the SMD force, in kJ mol-1

        :param path: Path to the file to which the data will be saved.
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
            np.save(outfile, self.smd_simulation_atom_positions)
            np.save(outfile, self.smd_simulation_work_done)

    def save_general_smd_data(self, path: str = None):
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

    def __init__(self, name: Optional[str] = None):

        super().__init__(name)

    def check_for_existing_smd_force(self):

        try:
            force_constant = self.simulation.context.getParameter("smd_k")
            n_forces = self.simulation.system.getNumForces()
            smd_force = self.simulation.system.getForce(n_forces - 1)
            params = smd_force.getParticleParameters(0)
            assert type(smd_force) == CustomExternalForce
            assert (
                smd_force.getEnergyFunction()
                == "0.5 * smd_k * periodicdistance(x, y, z, x0, y0, z0)^2"
                or smd_force.getEnergyFunction()
                == "0.5 * smd_k * ((x-x0)^2 + (y-y0)^2 + (z-z0)^2)"
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
        assert np.all(self.smd_simulation_atom_positions != 0.0)
        self._calculate_forces(self.smd_simulation_atom_positions)
        self._calculate_work_done()


class OpenMMSMDSimulationCOM(OpenMMSMDSimulation):
    """
    A class for performing constant velocity SMD on an OpenMM simulations,
    where the SMD force is applied to the centre of mass of a specified
    group of atoms.
    """

    def __init__(self, name: Optional[str] = None):

        super().__init__(name)
        self.com_positions: Optional[np.ndarray] = None

    def check_for_existing_smd_force(self):
        try:
            force_constant = self.simulation.context.getParameter("smd_k")
            n_forces = self.simulation.system.getNumForces()
            smd_force = self.simulation.system.getForce(n_forces - 1)
            assert type(smd_force) == CustomCentroidBondForce
            assert (
                smd_force.getEnergyFunction()
                == "0.5 * smd_k * pointdistance(x1, y1, z1, x0, y0, z0)^2"
            )
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

    def calculate_com(self, atom_positions: np.ndarray, atom_masses: np.ndarray):
        """
        Calculate the centre of mass of a group of N atoms.

        :param atom_positions: NumPy array of atom positions with dimensions (N, 3)
        :param atom_masses: NumPy array of atomic masses (AMU) with dimension (N)
        :return: NumPy array containing the position of the centre of mass of the atoms with dimension (3)
        """
        # TODO: Check this works and write tests!
        assert np.all(atom_positions.shape == np.array((self.smd_atom_indices.size, 3)))
        return np.sum(
            np.multiply(np.transpose(atom_positions), atom_masses), axis=1
        ) / np.sum(atom_masses)

    def calculate_com_trajectory(self):
        """
        Calculate the trajectory that the COM follows during the SMD simulation.
        """
        # TODO: Check this works and write tests!
        assert np.all(
            self.smd_simulation_atom_positions
            != np.zeros((self.smd_path.shape[0], self.smd_atom_indices.size, 3))
        )
        atom_masses = np.zeros(self.n_smd_atom_indices)
        for index in range(self.n_smd_atom_indices):
            atom_masses[index] = self.simulation.system.getParticleMass(
                self.smd_atom_indices[index]
            )._value
        self.com_positions = np.array(
            [
                self.calculate_com(self.smd_simulation_atom_positions[i], atom_masses)
                for i in range(self.smd_path.shape[0])
            ]
        )

    def calculate_cumulative_work_done(self):
        """
        Calculate the cumulative work done by the SMD force on the
        COM of the atoms with which it interacts over the SMD simulation.
        """
        self.calculate_com_trajectory()
        self._calculate_forces(self.com_positions)
        self._calculate_work_done()


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

    # TODO: Check that this works with PBCs correctly - may need to amend setUsesPBCs line...
    smd_force = CustomCentroidBondForce(
        1, "0.5 * smd_k * pointdistance(x1, y1, z1, x0, y0, z0)^2"
    )
    smd_force.addGlobalParameter("smd_k", force_constant)
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
        smd_force = CustomExternalForce(
            "0.5 * smd_k * periodicdistance(x, y, z, x0, y0, z0)^2"
        )
    else:
        smd_force = CustomExternalForce(
            "0.5 * smd_k * ((x-x0)^2 + (y-y0)^2 + (z-z0)^2)"
        )
    smd_force.addGlobalParameter("smd_k", force_constant)
    smd_force.addPerParticleParameter("x0")
    smd_force.addPerParticleParameter("y0")
    smd_force.addPerParticleParameter("z0")
    smd_force.setForceGroup(31)

    return smd_force
