"""
Minimum Narupa frame client that saves a trajectory to a file.

The client connects to a Narupa Frame server, and writes the frames it receives
into a file using MDAnalysis.

The client can be started either from python or from the command line. From the
command line, run this script with the `--help` option to see the usage. From
python, here is an example:

    from trajectory_client import write_trajectory_from_server
    write_trajectory_from_server('foo.xtc', address='localhost', port=54321)

This will connect to ``localhost:54321`` and writes the trajectory in 'foo.xtc'
in Gromacs's XTC format.

The client handles the trajectory formats that can be written by MDAnalysis;
the full list is available on
<https://www.mdanalysis.org/docs/documentation_pages/coordinates/init.html#supported-coordinate-formats>.
"""
import argparse

import numpy as np
import MDAnalysis as mda

from narupa.trajectory import FrameClient
from narupa.trajectory.frame_data import FrameData


def write_trajectory_from_server(destination, *, address: str, port: int):
    """
    Connect to a Narupa frame server and write the received frames into a file.

    :param destination: Path to the file to write. The format is guessed by
        MDAnalysis from the file extension.
    :param address: Host name to connect to.
    :param port: Port to connect to on the host.
    """
    with FrameClient(address=address, port=port) as client:
        print("Running...")
        frame_iter = client.subscribe_frames_iterate()
        frame_index, first_frame = next(frame_iter)
        first_positions = frame_to_ndarray(first_frame)
        universe = universe_from_positions(first_positions)
        writer_class = mda.coordinates.core.get_writer_for(destination,
                                                           multiframe=True)
        with writer_class(destination, n_atoms=first_positions.shape[0]) as writer:
            writer.write(universe.atoms)
            try:
                for frame_index, frame in frame_iter:
                    positions = frame_to_ndarray(frame)
                    universe.atoms.positions = positions
                    writer.write(universe.atoms)
            except KeyboardInterrupt:
                print(f"Finishing on frame {frame_index} due to keyboard interrupt.")


def universe_from_positions(positions: np.ndarray) -> mda.Universe:
    universe = mda.Universe.empty(n_atoms=positions.shape[0], trajectory=True)
    universe.atoms.positions = positions
    return universe


def frame_to_ndarray(frame: FrameData) -> np.ndarray:
    """
    Convert a frame received from the Narupa server to an array of coordinates.

    .. note::

        The frame is assumed to describes particles in a 3D space.

    :param frame: A frame received from a Narupa server.
    :return: A numpy array of ``np.float32`` with one row per atom, and 3
        columns. The coordinates are expressed in ångstöm.
    """
    # Here we assume that the coordinates are in 3 dimensions.
    # The frames obtained from Narupa use distances in nm, while MDAnalysis
    # expresses them in Å; this is where we do the unit conversion.
    positions = np.array(frame.particle_positions, dtype=np.float32) * 10
    return positions


def handle_user_args() -> argparse.Namespace:
    """
    Parse the arguments from the command line.

    :return: The namespace of arguments read from the command line.
    """
    description = "Connect to a Narupa trajectory server and write the trajectory."
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('--host', default=None)
    parser.add_argument('--port', type=int, default=None)
    parser.add_argument('destination')
    arguments = parser.parse_args()
    return arguments


def main():
    """
    Entry point for the command line.
    """
    arguments = handle_user_args()
    write_trajectory_from_server(
        arguments.destination,
        address=arguments.host,
        port=arguments.port,
    )


if __name__ == '__main__':
    main()
