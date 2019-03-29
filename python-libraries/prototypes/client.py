import grpc
from narupa.protocol.instance.molecule_provider_pb2_grpc import MoleculeProviderStub
from narupa.protocol.instance.get_frame_pb2 import GetFrameRequest, GetFrameResponse

import numpy as np
import MDAnalysis as mda

import argparse


class DummyClient:
    def __init__(self, destination, *, address: str, port: int):
        host = '{}:{}'.format(address, port)
        self.channel = grpc.insecure_channel(host)
        self.stub = MoleculeProviderStub(self.channel)
        self.destination = destination

    def write_trajectory(self):
        frame_iter = self.stub.SubscribeFrame(GetFrameRequest())
        first_frame = frame_to_ndarray(next(frame_iter))
        universe = mda.Universe.empty(n_atoms=first_frame.shape[0], trajectory=True)
        universe.atoms.positions = first_frame
        writer_class = mda.coordinates.core.get_writer_for(self.destination, multiframe=True)
        with writer_class(self.destination, n_atoms=first_frame.shape[0]) as writer:
            for frame in frame_iter:
                positions = frame_to_ndarray(frame)
                universe.atoms.positions = positions
                writer.write(universe.atoms)


def frame_to_ndarray(frame):
    raw_positions = (
        frame
        .frame.arrays['atom.position']
        .float_values
        .values
    )
    positions = np.array(raw_positions, dtype=np.float32).reshape((-1, 3)) * 10
    return positions


def handle_user_args():
    description = "Connect to a Narupa trajectory server and write the trajectory."
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('--host', default='localhost')
    parser.add_argument('--port', type=int, default=8000)
    parser.add_argument('destination')
    arguments = parser.parse_args()
    return arguments


def main():
    arguments = handle_user_args()
    dummy = DummyClient(
        arguments.destination,
        address=arguments.host,
        port=arguments.port,
    )
    dummy.write_trajectory()


if __name__ == '__main__':
    main()
