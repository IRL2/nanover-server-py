# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.
import time
from contextlib import contextmanager

import pytest
from ase.calculators.lj import LennardJones
from ase.md import VelocityVerlet

from nanover.ase.imd import NanoverASEDynamics
from nanover.ase.imd_calculator import ImdCalculator
from nanover.imd import ImdClient
from nanover.imd.particle_interaction import ParticleInteraction
from util import co_atoms, client_interaction
import numpy as np


@pytest.fixture
def interact_c():
    interaction = ParticleInteraction(
        position=[0, 1, 0],
        particles=[0],
        scale=20000.0,
        interaction_type="spring",
    )
    return interaction


@pytest.fixture
def interact_both():
    interaction = ParticleInteraction(
        position=[0, 1, 0],
        particles=[0, 1],
        scale=20000.0,
        interaction_type="spring",
        max_force=np.inf,
    )
    return interaction


@pytest.fixture
def imd_server_atoms_client():
    atoms = co_atoms()
    calculator = LennardJones()
    atoms.calc = calculator
    dynamics = VelocityVerlet(atoms, timestep=0.5)
    with NanoverASEDynamics.basic_imd(dynamics, port=0) as server:
        with ImdClient.insecure_channel(port=server.port) as client:
            yield server, atoms, client


def test_ase_imd_dynamics(imd_server_atoms_client):
    dynamics, atoms, imd_client = imd_server_atoms_client
    dynamics.run(5)


def test_ase_imd_dynamics_interaction(imd_server_atoms_client, interact_c):
    """
    Checks that an interactive force is integrated into the equations of motion correctly, by applying a huge
    interactive force and ensuring the resulting momentum is large in the direction of interaction.
    """
    dynamics, atoms, imd_client = imd_server_atoms_client

    imd_calculator = atoms.calc
    assert isinstance(imd_calculator, ImdCalculator)

    with client_interaction(imd_client, interact_c):
        time.sleep(0.1)
        assert len(imd_calculator.interactions) == 1
        dynamics.run(2)
        atom = dynamics.atoms[interact_c.particles[0]]
        assert atom.momentum[1] > 100


def test_ase_imd_dynamics_interaction_com(imd_server_atoms_client, interact_both):
    """
    Checks that an interactive force is integrated into the equations of motion correctly when applying a
    force to both atoms in the CO test system.
    """
    dynamics, atoms, imd_client = imd_server_atoms_client

    imd_calculator = atoms.calc
    assert isinstance(imd_calculator, ImdCalculator)

    with client_interaction(imd_client, interact_both):
        time.sleep(0.1)
        assert len(imd_calculator.interactions) == 1
        dynamics.run(1)
        for atom in dynamics.atoms:
            assert atom.momentum[1] > 50


def test_ase_imd_run_forever(imd_server_atoms_client):
    runner, atoms, imd_client = imd_server_atoms_client
    runner.run()
    time.sleep(0.1)
    runner.cancel_run(wait=True)
    number_of_steps = runner.dynamics.get_number_of_steps()
    assert number_of_steps > 0
    time.sleep(0.1)
    # check it really has stopped.
    assert number_of_steps == runner.dynamics.get_number_of_steps()


def test_get_calculator(imd_server_atoms_client):
    runner, atoms, imd_client = imd_server_atoms_client
    assert isinstance(runner.internal_calculator, LennardJones)
