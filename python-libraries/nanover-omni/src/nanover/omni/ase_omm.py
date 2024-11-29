import warnings
from typing import Optional

from ase import units, Atoms
from ase.md import Langevin
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from openmm.app import Simulation

from nanover.ase.omm_calculator import OpenMMCalculator
from nanover.omni.ase import ASESimulation


CONSTRAINTS_UNSUPPORTED_MESSAGE = (
    "The simulation contains constraints which will be ignored by this runner!"
)


class ASEOpenMMSimulation(ASESimulation):
    """
    A wrapper for ASE OpenMM simulations so they can be run inside the OmniRunner with some preset default dynamics.
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

    def __init__(self, name: Optional[str] = None):
        name = name or "Unnamed ASE OpenMM Simulation"

        super().__init__(name or "Unnamed ASE OpenMM Simulation")

        self.platform: Optional[str] = None
        self.simulation: Optional[Simulation] = None
        self.openmm_calculator: Optional[OpenMMCalculator] = None

    def load(self):
        """
        Load and set up the simulation if it isn't done already.
        """
        assert self.simulation is not None

        self.openmm_calculator = OpenMMCalculator(self.simulation)
        self.ase_atoms_to_frame_data = self.openmm_calculator.make_frame_converter()
        atoms = self.openmm_calculator.generate_atoms()

        self.initial_calc = self.openmm_calculator

        # we don't read this from the openmm xml
        self.dynamics = make_default_ase_omm_dynamics(atoms)

        self.atoms.calc = self.openmm_calculator

        # Set the momenta corresponding to T=300K
        MaxwellBoltzmannDistribution(self.atoms, temperature_K=300)

        if self.simulation.system.getNumConstraints() > 0:
            warnings.warn(CONSTRAINTS_UNSUPPORTED_MESSAGE)

        super().load()


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
