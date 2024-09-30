import pytest
from ase import units, Atoms
from ase.calculators.lj import LennardJones
from ase.md import VelocityVerlet

from nanover.omni.ase import ASESimulation

from common import app_server


@pytest.fixture
def example_ase(app_server, example_dynamics):
    sim = ASESimulation.from_dynamics(example_dynamics)
    sim.load()
    sim.reset(app_server)
    yield sim


@pytest.fixture
def example_dynamics():
    d = 1.1
    atoms = Atoms("CO", positions=[(0, 0, 0), (0, 0, d)], cell=[2, 2, 2], pbc=[1, 1, 1])
    calculator = LennardJones()
    atoms.calc = calculator
    dynamics = VelocityVerlet(atoms, timestep=0.5)
    yield dynamics


def test_step_interval(example_ase):
    """
    Test that advancing by one step increments the dynamics steps by frame_interval.
    """
    for i in range(5):
        assert (
            example_ase.dynamics.get_number_of_steps()
            == i * example_ase.frame_interval
        )
        example_ase.advance_by_one_step()
