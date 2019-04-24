"""Demonstrates IMD with ASE"""

from ase import units
from ase.calculators.emt import EMT
from ase.lattice.cubic import FaceCenteredCubic
from ase.md import Langevin
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution

from narupa.ase import NarupaASE
from narupa.ase.imd_calculator import ImdCalculator
from narupa.imd.imd_server import ImdServer
from narupa.trajectory import FrameServer


frame_server = FrameServer(address='localhost', port=54321)
imd_server = ImdServer(address='localhost', port=54322)

size = 2

# Set up a crystal
atoms = FaceCenteredCubic(directions=[[1, 0, 0], [0, 1, 0], [0, 0, 1]],
                          symbol="Cu",
                          size=(size, size, size),
                          pbc=True)

main_calculator = EMT()
imd_calculator = ImdCalculator(imd_server.service, main_calculator)

atoms.set_calculator(imd_calculator)

# Set the momenta corresponding to T=300K
MaxwellBoltzmannDistribution(atoms, 300 * units.kB)

# Room temperature simulation
dyn = Langevin(atoms, 1 * units.fs, units.kB * 300, 0.002)

# Attach the frame server.
dyn.attach(NarupaASE(atoms, frame_server), interval=1)
print("Starting dynamics")
while True:
    dyn.run(200)
