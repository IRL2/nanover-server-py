# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.
import time

import pytest
from ase.calculators.lj import LennardJones
from ase.md import VelocityVerlet

from narupa.ase.imd_server import ASEImdServer
from narupa.ase.imd_calculator import ImdCalculator
from narupa.core.timing import delayed_generator
from narupa.imd.particle_interaction import ParticleInteraction
from .util import co_atoms, imd_client


@pytest.fixture
def interact_c():
    interaction = ParticleInteraction(position=[0, 1, 0], particles=[0], scale=20000., interaction_type='spring')
    return interaction


@pytest.fixture
def interact_both():
    interaction = ParticleInteraction(position=[0, 1, 0], particles=[0, 1], scale=20000., interaction_type='spring')
    return interaction


@pytest.fixture
def imd():
    atoms = co_atoms()
    calculator = LennardJones()
    atoms.set_calculator(calculator)
    dynamics = VelocityVerlet(atoms, timestep=0.5)
    with ASEImdServer(dynamics) as imd:
        yield imd, atoms


def test_ase_imd_dynamics(imd):
    dynamics, atoms = imd
    dynamics.run(5)


def test_ase_imd_dynamics_interaction(imd, interact_c, imd_client):
    """
    Checks that an interactive force is integrated into the equations of motion correctly, by applying a huge
    interactive force and ensuring the resulting momentum is large in the direction of interaction.
    """
    dynamics, atoms = imd

    imd_calculator = atoms.get_calculator()
    assert isinstance(imd_calculator, ImdCalculator)

    imd_client.publish_interactions_async(delayed_generator([interact_c] * 10000, delay=0.025))
    time.sleep(0.1)
    assert len(imd_calculator.interactions) == 1

    dynamics.run(2)
    atom = dynamics.atoms[interact_c.particles[0]]
    assert atom.momentum[1] > 100


def test_ase_imd_dynamics_interaction_com(imd, interact_both, imd_client):
    """
    Checks that an interactive force is integrated into the equations of motion correctly when applying a
    force to both atoms in the CO test system.
    """
    dynamics, atoms = imd

    imd_calculator = atoms.get_calculator()
    assert isinstance(imd_calculator, ImdCalculator)
    imd_client.publish_interactions_async(delayed_generator([interact_both] * 10000, delay=0.025))
    time.sleep(0.1)
    assert len(imd_calculator.interactions) == 1

    dynamics.run(1)
    for atom in dynamics.atoms:
        assert atom.momentum[1] > 50


def test_ase_imd_run_forever(imd):
    runner, atoms = imd
    runner.run()
    time.sleep(0.1)
    runner.cancel_run(wait=True)
    number_of_steps = runner.dynamics.get_number_of_steps()
    assert number_of_steps > 0
    time.sleep(0.1)
    # check it really has stopped.
    assert number_of_steps == runner.dynamics.get_number_of_steps()


def test_get_calculator(imd):
    runner, atoms = imd
    assert isinstance(runner.internal_calculator, LennardJones)
