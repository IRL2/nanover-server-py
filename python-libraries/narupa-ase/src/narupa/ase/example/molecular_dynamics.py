"""Demonstrates molecular dynamics with constant energy."""

from ase import units
from ase.calculators.emt import EMT
from ase.lattice.cubic import FaceCenteredCubic
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.md.verlet import VelocityVerlet

from narupa.ase import NarupaASE
from narupa.trajectory import FrameServer

size = 2

# Set up a crystal
atoms = FaceCenteredCubic(directions=[[1, 0, 0], [0, 1, 0], [0, 0, 1]],
                          symbol="Cu",
                          size=(size, size, size),
                          pbc=True)

# Describe the interatomic interactions with the Effective Medium Theory
atoms.set_calculator(EMT())

# Set the momenta corresponding to T=300K
MaxwellBoltzmannDistribution(atoms, 300 * units.kB)

# We want to run MD with constant energy using the VelocityVerlet algorithm.
dyn = VelocityVerlet(atoms, 5 * units.fs)  # 5 fs time step.

server = FrameServer(address='localhost', port=54321)

# Now run the dynamics
dyn.attach(NarupaASE(atoms, server), interval=1)
print("Starting dynamics")
while True:
    dyn.run(200)
