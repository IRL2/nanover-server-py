# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.
"""
Helpers around useful potentials.
"""

from typing import Iterable
from simtk import openmm as mm
from simtk.unit import Quantity, kilojoule_per_mole, nanometer


def restraint_force(force_constant: Quantity = 100 * kilojoule_per_mole / nanometer ** 2):
    force = mm.CustomExternalForce('k * periodicdistance(x, y, z, x0, y0, z0)^2')
    force.addGlobalParameter('k', force_constant)
    force.addPerParticleParameter('x0')
    force.addPerParticleParameter('y0')
    force.addPerParticleParameter('z0')
    return force


def restrain_particles(
        force_constant: Quantity,
        positions: Quantity,
        particle_indices: Iterable[int],
):
    force = restraint_force(force_constant)
    for index in particle_indices:
        x0, y0, z0 = positions[index]
        force.addParticle(index, [x0, y0, z0])
    return force
