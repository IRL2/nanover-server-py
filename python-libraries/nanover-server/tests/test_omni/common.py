from contextlib import contextmanager
from pathlib import Path
from nanover.app import (
    NanoverImdApplication,
)
from nanover.app.types import AppServer
from nanover.omni import OmniRunner
from nanover.websocket import NanoverImdClient
from nanover.trajectory import FrameData

EXAMPLES_PATH = Path(__file__).parent
RECORDING_PATH = EXAMPLES_PATH / "nanotube-example-recording.nanover.zip"
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
    with NanoverImdClient.from_runner(runner) as client:
        yield client


@contextmanager
def make_connected_client_from_app_server(app_server: AppServer):
    with NanoverImdClient.from_app_server(app_server) as client:
        yield client


def connect_and_retrieve_first_frame_from_app_server(
    app_server: AppServer,
) -> FrameData:
    with make_connected_client_from_app_server(app_server) as client:
        return client.wait_until_first_frame()


@contextmanager
def make_app_server():
    with NanoverImdApplication.basic_server(port=0) as app_server:
        yield app_server
