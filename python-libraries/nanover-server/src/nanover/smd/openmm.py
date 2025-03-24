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
        name: Optional[str] = None
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
        name: Optional[str] = None
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

        self.n_smd_atom_indices: Optional[int] = None

        self.smd_force: Optional[
            Union[CustomExternalForce, CustomCentroidBondForce]
        ] = None

        self.checkpoint: Optional[Any] = None

        self.current_smd_force_position: Optional[np.ndarray] = None
        self.current_smd_force_position_index: Optional[int] = None
        self.smd_simulation_atom_positions: Optional[np.ndarray] = None


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


    def get_smd_atom_positions(self):
        """
        Retrieve the positions of the atoms with which the SMD force is
        interacting, and add them to the array of positions to save.
        """
        positions = self.simulation.context.getState(getPositions=True).getPositions(asNumpy=True)
        self.smd_simulation_atom_positions[self.current_smd_force_position_index] = positions[self.smd_atom_indices]



class OpenMMSMDSimulationAtom(OpenMMSMDSimulation):

    def __init__(self, name: Optional[str] = None):

        super().__init__(name)
        self.current_smd_force_position = self.smd_path[0]
        self.smd_simulation_atom_positions = np.zeros((self.smd_path.shape[0], 3))


    def add_smd_force_to_system(self):

        x0, y0, z0 = self.current_smd_force_position
        smd_force = smd_single_atom_force(self.smd_force_constant)
        smd_force.addParticle(self.smd_atom_indices, [x0, y0, z0])
        self.smd_force = smd_force
        self.simulation.system.addForce(self.smd_force)


    def update_smd_force_position(self):

        x0, y0, z0 = self.current_smd_force_position
        self.smd_force.setParticleParameters(0, self.smd_atom_indices, [x0, y0, z0])
        self.smd_force.updateParametersInContext(self.simulation.context)



class OpenMMSMDSimulationCOM(OpenMMSMDSimulation):

    def __init__(self, name: Optional[str] = None):

        super().__init__(name)
        self.current_smd_force_position = self.smd_path[0]
        self.smd_simulation_atom_positions = np.zeros((self.smd_path.shape[0], self.n_smd_atom_indices, 3))


    def add_smd_force_to_system(self):

        x0, y0, z0 = self.current_smd_force_position
        smd_force = smd_com_force(self.smd_force_constant)
        smd_force.addGroup(self.smd_atom_indices)
        smd_force.addBond([0], [x0, y0, z0])
        self.smd_force = smd_force
        self.simulation.system.addForce(self.smd_force)


    def update_smd_force_position(self):

        x0, y0, z0 = self.current_smd_force_position
        self.smd_force.setBondParameters(0, 0, [x0, y0, z0])
        self.smd_force.updateParametersInContext(self.simulation.context)


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
