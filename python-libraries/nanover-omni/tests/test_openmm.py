import pytest

from nanover.omni.openmm import OpenMMSimulation

from common import app_server, ARGON_XML_PATH
from nanover.openmm.imd import add_imd_force_to_system


@pytest.fixture
def example_openmm(app_server):
    sim = OpenMMSimulation.from_xml_path(ARGON_XML_PATH)
    sim.load()
    sim.reset(app_server)
    yield sim


def test_simulation_without_imd_force(example_openmm):
    """
    Test that creation from a simulation without imd force fails.
    """
    simulation = example_openmm.simulation
    simulation.system.removeForce(simulation.system.getNumForces() - 1)

    with pytest.raises(ValueError):
        OpenMMSimulation.from_simulation(simulation)


def test_simulation_with_multiple_imd_force(example_openmm):
    """
    Test that creation from a simulation with multiple imd forces fails.
    """
    simulation = example_openmm.simulation

    # The forces added to the system will not be accounted for when running
    # the dynamics until the context is reset as the system is already
    # compiled in a context. It does not matter here, as the force is still
    # listed in the system, which is what we check.
    add_imd_force_to_system(simulation.system)

    with pytest.warns(UserWarning, match="force"):
        OpenMMSimulation.from_simulation(simulation)


def test_step_interval(example_openmm):
    """
    Test that advancing by one step increments the dynamics steps by frame_interval.
    """
    for i in range(5):
        assert (
            example_openmm.simulation.currentStep == i * example_openmm.frame_interval
        )
        example_openmm.advance_by_one_step()
