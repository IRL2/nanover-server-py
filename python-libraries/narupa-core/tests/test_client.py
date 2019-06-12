import grpc
import pytest

from .test_frame_server import simple_frame_data, frame_server
from .imd.test_imd_service import interaction
from narupa.app.client import NarupaClient


def test_no_imd(frame_server, interaction):
    client = NarupaClient(run_imd=False)
    with pytest.raises(grpc.RpcError):
        client.start_interaction(interaction)