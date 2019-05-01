import time

import pytest
from ase import Atoms
from ase.calculators.lj import LennardJones
from ase.md import Langevin

from narupa.ase.ase_imd import ASEImd
from narupa.ase.imd_calculator import ImdCalculator
from narupa.imd.imd_client import ImdClient, delayed_generator
from narupa.imd.imd_server import ImdServer
from narupa.imd.interaction import Interaction


def co_atoms():
    d = 1.1
    co = Atoms('CO', positions=[(0, 0, 0), (0, 0, d)],
               cell=[2, 2, 2],
               pbc=[1, 1, 1])
    return co


@pytest.fixture
def atoms():
    return co_atoms()


@pytest.fixture
def interact_c():
    interaction = Interaction(position=[1, 0, 0], particles=[0], scale=100., interaction_type='spring')
    return interaction


@pytest.fixture
def imd():
    atoms = co_atoms()
    calculator = LennardJones()
    atoms.set_calculator(calculator)
    dynamics = Langevin(atoms, timestep=0.5, temperature=300, friction=1.0)
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
    dynamics, atoms = imd

    imd_calculator = atoms.get_calculator()
    assert isinstance(imd_calculator, ImdCalculator)
    imd_client.publish_interactions_async(delayed_generator([interact_c] * 50, delay=0.025))
    time.sleep(0.05)
    assert len(imd_calculator.interactions) == 1

    dynamics.run(1)
    atom = dynamics.atoms[interact_c.particles[0]]
    print(atom)

