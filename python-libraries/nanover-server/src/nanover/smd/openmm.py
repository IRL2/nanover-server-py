from os import PathLike
from pathlib import Path
from typing import Optional, Union, Any
from abc import abstractmethod

import numpy as np

from openmm import CustomExternalForce, CustomCentroidBondForce
from openmm.app import Simulation

from nanover.openmm import serializer


class OpenMMSMDSimulation:
    """
    A wrapper for performing constant velocity SMD on an OpenMM simulation.
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
        """

        # Create instance of SMD simulation based on type of indices passed
        assert smd_atom_indices.size >= 1
        if smd_atom_indices.size > 1:
            sim = super(cls, OpenMMSMDSimulationCOM).__new__(OpenMMSMDSimulationCOM)
        else:
            sim = super(cls, OpenMMSMDSimulationAtom).__new__(OpenMMSMDSimulationAtom)

        sim.name = name
        sim.simulation = simulation
        sim.smd_atom_indices = smd_atom_indices
        sim.smd_path = smd_path
        sim.smd_force_constant = smd_force_constant
        sim.n_smd_atom_indices = sim.smd_atom_indices.size
        sim.define_smd_simulation_atom_positions_array()

        # Create SMD force and add it to the system
        sim.current_smd_force_position = sim.smd_path[0]
        sim.current_smd_force_position_index = 0
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
        sim.smd_atom_indices = smd_atom_indices
        sim.n_smd_atom_indices = sim.smd_atom_indices.size
        sim.smd_path = smd_path
        sim.smd_force_constant = smd_force_constant
        sim.define_smd_simulation_atom_positions_array()

        # Create SMD force and add it to the system
        sim.current_smd_force_position = sim.smd_path[0]
        sim.current_smd_force_position_index = 0
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

    @abstractmethod
    def define_smd_simulation_atom_positions_array(self):
        """
        Define the array to which the positions of the atoms with which the
        SMD force interacts will be saved over the course of the SMD simulation.
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
            raise AssertionError("Restraint is not located at the initial position "
                                 "of the SMD force, equilibration aborted. To prepare the system "
                                 "appropriately, the restraint should be placed at the first point "
                                 "on the path that the SMD force will take during the SMD simulation.")

        self.simulation.step(n_steps)

    def run_smd(self, progress_interval: Optional[int] = 100000):
        """
        Perform an SMD simulation on the system using the SMD force and
        path defined, store the positions of the atoms with which the
        SMD force interacts, and calculate the work done along the reaction
        coordinate.
        """

        # Get initial atom positions for SMD force at initial position
        assert(self.smd_simulation_atom_positions is not None and np.all(
            self.current_smd_force_position == self.smd_path[0]) and
            self.current_smd_force_position_index == 0
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


    def calculate_forces(self, interaction_centre_positions):
        """
        Calculate the SMD forces that acted on the system during the simulation.
        """
        assert np.all(self.smd_path.shape == interaction_centre_positions.shape)
        self.smd_simulation_forces = -self.smd_force_constant * (
            interaction_centre_positions - self.smd_path
        )


    def _calculate_work_done(self):
        """
        Calculate the cumulative work done along the reaction coordinate
        """
        assert self.smd_path is not None and self.smd_simulation_forces is not None
        smd_force_displacements = np.diff(self.smd_path, axis=0)
        work_done_array = np.zeros(self.smd_simulation_forces.shape[0])
        for i in range(smd_force_displacements.shape[0]):
            work_done_array[i+1] = np.dot(
                self.smd_simulation_forces[i], smd_force_displacements[i]
            )
        self.smd_simulation_work_done = np.cumsum(work_done_array, axis=0)


    def save_smd_simulation_data(self, path: str = None):

        if path is None:
            raise ValueError("Output file path cannot be None. Please specify an output file path.")

        elif self.smd_simulation_work_done is None:
            raise ValueError("Missing values for the work done. This data can only be saved after"
                             "the SMD calculation is completed.")

        elif self.smd_simulation_atom_positions is None or np.all(self.smd_simulation_atom_positions == 0.0):
            raise ValueError("Missing values for the atom positions. This data can only be saved after"
                             "the SMD calculation is completed.")

        with open(path, "wb") as outfile:
            np.save(outfile, self.smd_simulation_atom_positions)
            np.save(outfile, self.smd_simulation_work_done)


    def save_general_smd_data(self, path: str = None):

        if path is None:
            raise ValueError("Output file path cannot be None. Please specify an output file path.")

        with open(path, "wb") as outfile:
            np.save(outfile, self.smd_atom_indices)
            np.save(outfile, self.smd_path)
            np.save(outfile, self.smd_force_constant)
            np.save(outfile, self.simulation.integrator.getTemperature()._value)
            np.save(outfile, self.simulation.integrator.getStepSize()._value)


class OpenMMSMDSimulationAtom(OpenMMSMDSimulation):

    def __init__(self, name: Optional[str] = None):

        super().__init__(name)

    def add_smd_force_to_system(self):

        x0, y0, z0 = self.smd_path[0]
        smd_force = smd_single_atom_force(self.smd_force_constant)
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
        self.calculate_forces(self.smd_simulation_atom_positions)
        self._calculate_work_done()


class OpenMMSMDSimulationCOM(OpenMMSMDSimulation):

    def __init__(self, name: Optional[str] = None):

        super().__init__(name)
        self.com_positions: Optional[np.ndarray] = None

    def add_smd_force_to_system(self):

        x0, y0, z0 = self.smd_path[0]
        smd_force = smd_com_force(self.smd_force_constant)
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
        Calculate the centre of mass of a group of atoms.
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
        self.calculate_forces(self.com_positions)
        self._calculate_work_done()


def smd_com_force(force_constant: float):
    """
    Defines a harmonic restraint force for the COM of a group of atoms for performing SMD.
    """

    # TODO: Check that this works with PBCs correctly - may need to amend setUsesPBCs line...
    smd_force = CustomCentroidBondForce(
        1, "0.5 * smd_k * pointdistance(x1, y1, z1, x0, y0, z0)^2"
    )
    smd_force.addGlobalParameter("smd_k", force_constant)
    smd_force.addPerBondParameter("x0")
    smd_force.addPerBondParameter("y0")
    smd_force.addPerBondParameter("z0")
    smd_force.setUsesPeriodicBoundaryConditions(True)
    smd_force.setForceGroup(31)

    return smd_force


def smd_single_atom_force(force_constant: float):
    """
    Defines a harmonic restraint force for a single atom for performing SMD.
    """

    smd_force = CustomExternalForce(
        "0.5 * smd_k * periodicdistance(x, y, z, x0, y0, z0)^2"
    )
    smd_force.addGlobalParameter("smd_k", force_constant)
    smd_force.addPerParticleParameter("x0")
    smd_force.addPerParticleParameter("y0")
    smd_force.addPerParticleParameter("z0")
    smd_force.setForceGroup(31)

    return smd_force
