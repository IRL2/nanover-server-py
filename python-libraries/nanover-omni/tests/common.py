from contextlib import contextmanager
from pathlib import Path
from nanover.app import (
    NanoverImdClient,
    NanoverApplicationServer,
    NanoverImdApplication,
)
from nanover.omni import OmniRunner

EXAMPLES_PATH = Path(__file__).parent
RECORDING_PATH_TRAJ = EXAMPLES_PATH / "nanotube-example-recording.traj"
RECORDING_PATH_STATE = EXAMPLES_PATH / "nanotube-example-recording.state"
ARGON_XML_PATH = EXAMPLES_PATH / "argon_simulation.xml"


@contextmanager
def make_loaded_sim(sim):
    with make_app_server() as app_server:
        sim.load()
        sim.reset(app_server)
        yield sim


@contextmanager
def make_loaded_sim_with_interactions(sim, *interactions):
    with make_app_server() as app_server:
        sim.load()
        sim.reset(app_server)

        for i, interaction in enumerate(interactions):
            app_server.imd.insert_interaction(f"interaction.{i}", interaction)

        yield app_server, sim


@contextmanager
def make_runner(*simulations):
    with OmniRunner.with_basic_server(*simulations, port=0) as runner:
        yield runner


@contextmanager
def make_connected_client_from_runner(runner):
    with make_connected_client_from_app_server(runner.app_server) as client:
        yield client


@contextmanager
def make_connected_client_from_app_server(app_server: NanoverApplicationServer):
    with NanoverImdClient.connect_to_single_server(port=app_server.port) as client:
        yield client


def connect_and_retrieve_first_frame_from_app_server(
    app_server: NanoverApplicationServer,
):
    with make_connected_client_from_app_server(app_server) as client:
        client.subscribe_to_frames()
        return client.wait_until_first_frame()


@contextmanager
def make_app_server():
    with NanoverImdApplication.basic_server(port=0) as app_server:
        yield app_server
