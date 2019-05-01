import time

import pytest
from ase import Atoms
from ase.calculators.lj import LennardJones
from ase.md import Langevin, VelocityVerlet

from narupa.ase.ase_imd import ASEImd
from narupa.ase.imd_calculator import ImdCalculator
from narupa.imd.imd_client import ImdClient, delayed_generator
from narupa.imd.imd_server import ImdServer
from narupa.imd.interaction import Interaction


def co_atoms():
    d = 1.1
    co = Atoms('CO', positions=[(0, 0, 0), (0, 0, d)],
               cell=[20, 20, 20],
               pbc=[1, 1, 1])
    return co


@pytest.fixture
def atoms():
    return co_atoms()


@pytest.fixture
def interact_c():
    interaction = Interaction(position=[0, 1, 0], particles=[0], scale=20000., interaction_type='spring')
    return interaction

@pytest.fixture
def interact_both():
    interaction = Interaction(position=[0, 1, 0], particles=[0,1], scale=20000., interaction_type='spring')
    return interaction


@pytest.fixture
def imd():
    atoms = co_atoms()
    calculator = LennardJones()
    atoms.set_calculator(calculator)
    dynamics = VelocityVerlet(atoms, timestep=0.5)
    imd = ASEImd(dynamics)
    yield imd, atoms
    imd.close()

@pytest.fixture
def imd_client():
    client = ImdClient(address='localhost', port=54322)
    yield client
    client.close()

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

    dynamics.run(10)
    atom = dynamics.atoms[interact_c.particles[0]]
    assert atom.momentum[1] > 200

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

    dynamics.run(10)
    for atom in dynamics.atoms:
        assert atom.momentum[1] > 200

