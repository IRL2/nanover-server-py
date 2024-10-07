import pytest
from ase import units, Atoms
from ase.calculators.lj import LennardJones
from ase.md import VelocityVerlet

from nanover.imd import ParticleInteraction
from nanover.omni import OmniRunner
from nanover.omni.ase import ASESimulation

from common import app_server


@pytest.fixture
def example_ase(app_server, example_dynamics):
    sim = ASESimulation.from_ase_dynamics(example_dynamics)
    sim.load()
    sim.reset(app_server)
    yield sim


@pytest.fixture
def example_dynamics():
    atoms = Atoms("C", positions=[(0, 0, 0)], cell=[2, 2, 2])
    atoms.calc = LennardJones()
    dynamics = VelocityVerlet(atoms, timestep=0.5)
    yield dynamics


def test_step_interval(example_ase):
    """
    Test that advancing by one step increments the dynamics steps by frame_interval.
    """
    for i in range(5):
        assert (
            example_ase.dynamics.get_number_of_steps() == i * example_ase.frame_interval
        )
        example_ase.advance_by_one_step()


def test_dynamics_interaction(example_ase):
    """
    Test that example dynamics responds to interactions.
    """

    with OmniRunner.with_basic_server(example_ase, port=0) as runner:
        runner.app_server.imd.insert_interaction(
            "interaction.0",
            ParticleInteraction(
                position=(0.0, 0.0, 10.0),
                particles=[0],
                interaction_type="constant",
            ),
        )
        runner.next()
        runner.pause()
        for _ in range(30):
            example_ase.advance_by_one_step()

    positions = example_ase.atoms.get_positions()
    (x, y, z) = positions[0]

    assert z > 2.5
