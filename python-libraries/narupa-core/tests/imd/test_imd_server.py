from typing import Generator, Tuple

import grpc
import pytest
from narupa.imd.imd_client import ImdClient
from narupa.imd.imd_server import ImdServer
from narupa.imd.particle_interaction import ParticleInteraction
from narupa.protocol.imd.imd_pb2_grpc import InteractiveMolecularDynamicsStub


@pytest.fixture
def imd_server() -> Generator[ImdServer, None, None]:
    with ImdServer(address='localhost', port=0) as server:
        yield server


@pytest.fixture
def imd_server_client(imd_server) -> Generator[Tuple[ImdServer, ImdClient], None, None]:
    with ImdClient.insecure_channel(address='localhost', port=imd_server.port) as client:
        yield imd_server, client


@pytest.fixture
def imd_server_stub(imd_server) -> Generator[Tuple[ImdServer, InteractiveMolecularDynamicsStub], None, None]:
    channel = grpc.insecure_channel(f"localhost:{imd_server.port}")
    try:
        stub = InteractiveMolecularDynamicsStub(channel)
        yield imd_server, stub
    finally:
        channel.close()


@pytest.fixture
def interaction():
    return ParticleInteraction()



