from os import PathLike
from pathlib import Path
from typing import Optional, Union, Any

import numpy as np

from openmm import CustomExternalForce, CustomCentroidBondForce
from openmm.app import Simulation

from nanover.openmm import serializer

class OpenMMSMDSimulation:
    """
    A wrapper for performing constant velocity SMD on an OpenMM simulation.
    """

    @classmethod
    def from_simulation(cls, simulation: Simulation, smd_atom_indices: np.ndarray, smd_path: np.ndarray, smd_force_constant: float, *, name: Optional[str] = None):
        """
        Construct the SMD simulation from an existing OpenMM simulation.
        """
        sim = cls(name)
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
    def from_xml_path(cls, path: PathLike[str], smd_atom_indices: np.ndarray, smd_path: np.ndarray, smd_force_constant: float, *, name: Optional[str] = None):
        """
        Construct the SMD simulation from an existing NanoVer OpenMM XML file located at a given path.
        """
        if name is None:
            name = Path(path).stem
        sim = cls(name)
        sim.xml_path = path

        # Load the simulation from the path
        with open(sim.xml_path) as infile:
            sim.simulation = serializer.deserialize_simulation(infile)

        sim.smd_atom_indices = smd_atom_indices
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

        self.smd_force: Optional[Union[CustomExternalForce, CustomCentroidBondForce]] = None

        self.checkpoint: Optional[Any] = None


    def add_smd_force_to_system(self):
        """
        Add the appropriate SMD force to the system.
        """
        assert(self.simulation is not None
               and self.smd_force_constant is not None
               and self.smd_path is not None
               and self.smd_atom_indices is not None
        )

        self.smd_force = define_smd_force(force_constant=self.smd_force_constant,
                                     atom_indices=self.smd_atom_indices,
                                     restraint_position=self.smd_path[0])

        self.simulation.system.addForce(self.smd_force)


def smd_com_force(force_constant: float):
    """
    Defines a harmonic restraint force for the COM of a group of atoms for performing SMD.
    """

    # TODO: Check that this works with PBCs correctly - may need to amend setUsesPBCs line...
    smd_force = CustomCentroidBondForce(1, "0.5 * smd_k * pointdistance(x, y, z, x0, y0, z0)^2")
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

    smd_force = CustomExternalForce("0.5 * smd_k * periodicdistance(x, y, z, x0, y0, z0)^2")
    smd_force.addGlobalParameter("smd_k", force_constant)
    smd_force.addPerParticleParameter("x0")
    smd_force.addPerParticleParameter("y0")
    smd_force.addPerParticleParameter("z0")
    smd_force.setForceGroup(31)

    return smd_force


def define_smd_force(force_constant: float, atom_indices: np.ndarray, restraint_position: np.ndarray):
    """
    Define a harmonic restraint force for performing SMD.
    """

    x0, y0, z0 = restraint_position

    try:
        # If multiple atom indices are passed, create an SMD force for the COM of the atoms
        assert(atom_indices.shape[0] > 1)
        smd_force = smd_com_force(force_constant)
        smd_force.addGroup(atom_indices)
        smd_force.addBond([0], [x0, y0, z0])

    except IndexError:
        # If only single index is passed, create SMD force for single atom
        smd_force = smd_single_atom_force(force_constant)
        smd_force.addParticle(atom_indices, [x0, y0, z0])

    return smd_force
