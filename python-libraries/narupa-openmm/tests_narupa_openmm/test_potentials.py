# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.

from simtk import openmm as mm
from narupa.openmm import potentials
from .simulation_utils import basic_simulation
from simtk.unit import kilojoule_per_mole, nanometer


def test_restraint_force_type():
    force = potentials.restraint_force()
    assert isinstance(force, mm.CustomExternalForce)


def test_restraint_force_global_parameters():
    force = potentials.restraint_force()
    global_parameters = [
        force.getGlobalParameterName(i)
        for i in range(force.getNumGlobalParameters())
    ]
    assert global_parameters == ['k']


def test_restraint_force_per_particle_parameters():
    force = potentials.restraint_force()
    per_particle_parameters = [
        force.getPerParticleParameterName(i)
        for i in range(force.getNumPerParticleParameters())
    ]
    assert per_particle_parameters == ['x0', 'y0', 'z0']


def test_restrain_particles(basic_simulation):
    positions = basic_simulation.context.getState(getPositions=True).getPositions()
    force = potentials.restrain_particles(
        positions=positions,
        particle_indices=[0, 4],
    )

    per_particle_parameters = [
        force.getParticleParameters(i)
        for i in range(force.getNumParticles())
    ]
    expectation = [
        [0, tuple(positions[0].value_in_unit(nanometer))],
        [4, tuple(positions[4].value_in_unit(nanometer))],
    ]

    assert force.getNumParticles() == 2
    assert per_particle_parameters == expectation

