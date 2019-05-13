# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.
"""
Demonstrates interactive molecular dynamics for a small EMT crystal.
"""

from ase import units
from ase.calculators.emt import EMT
from ase.lattice.cubic import FaceCenteredCubic
from ase.md import Langevin
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution

from narupa.ase.imd_server import ASEImdServer

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

dyn = Langevin(atoms, 1 * units.fs, 300, 0.1)

imd = ASEImdServer(dyn)
while True:
    imd.run(100)
