import warnings
from os import PathLike
from pathlib import Path
from typing import Optional

from ase import units, Atoms
from ase.md import Langevin, MDLogger
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from openmm.app import Simulation

from nanover.app import NanoverImdApplication
from nanover.ase.imd_calculator import ImdCalculator
from nanover.ase.openmm import OpenMMCalculator
from nanover.ase.openmm.frame_adaptor import (
    openmm_ase_atoms_to_frame_data,
)
from nanover.omni.ase import ASESimulation
from nanover.openmm import serializer


CONSTRAINTS_UNSUPPORTED_MESSAGE = (
    "The simulation contains constraints which will be ignored by this runner!"
)


class ASEOpenMMSimulation(ASESimulation):
    """
    A wrapper for ASE OpenMM simulations so they can be run inside the OmniRunner.
    """

    @classmethod
    def from_simulation(
        cls,
        simulation: Simulation,
        *,
        name: Optional[str] = None,
    ):
        """
        Construct this from an existing ASE OpenMM simulation.
        :param simulation: An existing ASE OpenMM Simulation
        :param name: An optional name for the simulation instead of default
        """
        sim = cls(name)
        sim.simulation = simulation
        return sim

    @classmethod
    def from_xml_path(cls, path: PathLike[str], *, name: Optional[str] = None):
        """
        Construct this from an existing NanoVer OpenMM XML file at a given path.
        :param path: Path of the NanoVer OpenMM XML file
        :param name: An optional name for the simulation instead of filename
        """
        sim = cls(name or Path(path).stem)
        sim.xml_path = path
        return sim

    def __init__(self, name: Optional[str] = None):
        name = name or "Unnamed ASE OpenMM Simulation"

        super().__init__(name or "Unnamed ASE OpenMM Simulation")

        self.ase_atoms_to_frame_data = openmm_ase_atoms_to_frame_data

        self.xml_path: Optional[PathLike[str]] = None

        self.platform: Optional[str] = None
        self.simulation: Optional[Simulation] = None
        self.openmm_calculator: Optional[OpenMMCalculator] = None

    def load(self):
        """
        Load and set up the simulation if it isn't done already.
        """
        if self.simulation is None and self.xml_path is not None:
            with open(self.xml_path) as infile:
                self.simulation = serializer.deserialize_simulation(
                    infile, platform_name=self.platform
                )

        assert self.simulation is not None

        self.openmm_calculator = OpenMMCalculator(self.simulation)
        atoms = self.openmm_calculator.generate_atoms()

        # we don't read this from the openmm xml
        self.dynamics = make_default_ase_omm_dynamics(atoms)

        self.atoms.calc = self.openmm_calculator

        # Set the momenta corresponding to T=300K
        MaxwellBoltzmannDistribution(self.atoms, temperature_K=300)

        super().load()

    def reset(self, app_server: NanoverImdApplication):
        """
        Reset the simulation to its initial conditions, reset IMD interactions, and reset frame stream to begin with
        topology and continue.
        :param app_server: The app server hosting the frame publisher and imd state
        """
        assert (
            self.simulation is not None
            and self.dynamics is not None
            and self.atoms is not None
            and self.checkpoint is not None
            and self.openmm_calculator is not None
        )

        self.app_server = app_server

        self.atoms.set_positions(self.checkpoint.positions)
        self.atoms.set_velocities(self.checkpoint.velocities)
        self.atoms.set_cell(self.checkpoint.cell)

        if self.simulation.system.getNumConstraints() > 0:
            warnings.warn(CONSTRAINTS_UNSUPPORTED_MESSAGE)

        self.atoms.calc = ImdCalculator(
            self.app_server.imd,
            self.openmm_calculator,
            dynamics=self.dynamics,
        )

        # send the initial topology frame
        frame_data = self.make_topology_frame()
        self.app_server.frame_publisher.send_frame(0, frame_data)
        self.frame_index = 1

        # TODO: deal with this when its clear if dynamics should be reconstructed or not..
        if self.verbose:
            self.dynamics.attach(
                MDLogger(
                    self.dynamics,
                    self.atoms,
                    "-",
                    header=True,
                    stress=False,
                    peratom=False,
                ),
                interval=100,
            )


def make_default_ase_omm_dynamics(atoms: Atoms):
    # We do not remove the center of mass (fixcm=False). If the center of
    # mass translations should be removed, then the removal should be added
    # to the OpenMM system.
    dynamics = Langevin(
        atoms=atoms,
        timestep=1 * units.fs,
        temperature_K=300,
        friction=1e-2,
        fixcm=False,
    )

    return dynamics
