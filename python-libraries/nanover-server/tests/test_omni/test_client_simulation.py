from nanover.app import OmniRunner
from nanover.testing import assert_equal_soon
from nanover.websocket.client.app_client import NanoverImdClient
from .test_openmm import make_example_openmm


def test_client_simulation_topology():
    """
    Test that an OmniRunner made from one client can be controlled via another and result in the expected topology
    seen on the controlling client.
    """
    simulations = [make_example_openmm()]

    with (
        OmniRunner.with_basic_server(port=0) as empty_runner,
        NanoverImdClient.from_runner(empty_runner) as simulating_client,
        NanoverImdClient.from_runner(empty_runner) as observing_client,
        OmniRunner.from_client(simulating_client, *simulations),
    ):
        for i, simulation in enumerate(simulations):
            observing_client.run_load(index=i)
            assert_equal_soon(
                lambda: observing_client.current_frame.particle_elements.all(),
                lambda: simulation.make_topology_frame().particle_elements.all(),
            )
