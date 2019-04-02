"""
Minimum Narupa frame client that saves a trajectory to a file.

The client connects to a Narupa Frame server, and writes the frames it receives
into a file using MDAnalysis.

The client can be started either from python or from the command line. From the
command line, run this script with the `--help` option to see the usage. From
python, here is an example::

    from client import DummyClient
    trajectory_client = DummyClient('foo.xtc', address='localhost', port=9000)
    trajectory_client.write_trajectory()

This will connect to ``localhost:9000`` and writes the trajectory in 'foo.xtc'
in Gromacs's XTC format.

The client handles the trajectory formats that can be written by MDAnalysis;
the full list is available on
<https://www.mdanalysis.org/docs/documentation_pages/coordinates/init.html#supported-coordinate-formats>.
"""
import argparse

import numpy as np
import MDAnalysis as mda

import grpc
from narupa.protocol.instance.molecule_provider_pb2_grpc import MoleculeProviderStub
from narupa.protocol.instance.get_frame_pb2 import GetFrameRequest, GetFrameResponse


# At this stage, this should not be a class. A class with a single public
# method should be a function. This is still a class because I over-think the
# future and assumes there will be other methods latter...
class DummyClient:
    """
    Connect to a Narupa frame server and write the received frames into a file.

    :param destination: Path to the file to write. The format is guessed by
        MDAnalysis from the file extension.
    :param address: Host name to connect to.
    :param port: Port to connect to on the host.
    """
    def __init__(self, destination, *, address: str, port: int):
        host = '{}:{}'.format(address, port)
        self.channel = grpc.insecure_channel(host)
        self.stub = MoleculeProviderStub(self.channel)
        self.destination = destination

    def write_trajectory(self):
        """
        Start receiving frames from the server and writing to file.
        """
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


def frame_to_ndarray(frame: GetFrameResponse) -> np.ndarray:
    """
    Convert a frame received from the Narupa server to an array of coordinates.

    .. note::

        The frame is assumed to describes particles in a 3D space.

    :param frame: A frame received from a Narupa server.
    :return: A numpy array of ``np.float32`` with one row per atom, and 3
        columns. The coordinates are expressed in ångstöm.
    """
    raw_positions = (
        frame
        .frame.arrays['atom.position']
        .float_values
        .values
    )
    # Here we assume that the coordinates are in 3 dimensions.
    # The frames obtained from Narupa use distances in nm, while MDAnalysis
    # expresses them in Å; this is where we do the unit conversion.
    positions = np.array(raw_positions, dtype=np.float32).reshape((-1, 3)) * 10
    return positions


def handle_user_args() -> argparse.Namespace:
    """
    Parse the arguments from the command line.

    :return: The namespace of arguments read from the command line.
    """
    description = "Connect to a Narupa trajectory server and write the trajectory."
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('--host', default='localhost')
    parser.add_argument('--port', type=int, default=8000)
    parser.add_argument('destination')
    arguments = parser.parse_args()
    return arguments


def main():
    """
    Entry point for the command line.
    """
    arguments = handle_user_args()
    dummy = DummyClient(
        arguments.destination,
        address=arguments.host,
        port=arguments.port,
    )
    dummy.write_trajectory()


if __name__ == '__main__':
    main()
