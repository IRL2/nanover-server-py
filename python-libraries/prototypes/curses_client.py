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

import math
import numpy as np
import MDAnalysis as mda
import curses

import grpc
from narupa.protocol.trajectory import (
    TrajectoryServiceStub,
    GetFrameRequest,
    GetFrameResponse,
)
from narupa.trajectory.frame_data import POSITIONS

cells = [" ", ".", "o", "O", "@"]

def render_frame(stdscr, positions):
    stdscr.clear()
    h, w = stdscr.getmaxyx()

    h -= 1
    w -= 1

    xmin = positions[:,0].min()
    ymin = positions[:,1].min()

    xmax = positions[:,0].max()
    ymax = positions[:,1].max()

    xr = xmax - xmin
    yr = ymax - ymin

    scale = min(w / xr, h / yr)

    ox = (w - xr * scale) / 2
    oy = (h - yr * scale) / 2

    positions[:,0] -= xmin
    positions[:,0] *= scale
    positions[:,0] += ox

    positions[:,1] -= ymin
    positions[:,1] *= scale
    positions[:,1] -= oy

    counts = {}

    for position in positions:
        coord = int(round(position[0])), int(round(position[1]))
        count = counts[coord] if coord in counts else 0
        counts[coord] = count + 1

    for coord, count in counts.items():
        x, y = coord

        if x < 0 or x >= w or y < 0 or y >= h:
            continue

        cell = cells[min(math.ceil(count / 3), 4)]
        stdscr.addstr(y, x, cell)

    stdscr.refresh()

def write_trajectory_from_server(stdscr, *, address: str, port: int):
    """
    Connect to a Narupa frame server and write the received frames into a file.

    :param destination: Path to the file to write. The format is guessed by
        MDAnalysis from the file extension.
    :param address: Host name to connect to.
    :param port: Port to connect to on the host.
    """
    stdscr.clear()
    stdscr.nodelay(True)
    stdscr.addstr(0, 0, "Connecting...")
    stdscr.refresh()

    host = '{}:{}'.format(address, port)
    channel = grpc.insecure_channel(host)
    stub = TrajectoryServiceStub(channel)
    frame_iter = stub.SubscribeFrames(GetFrameRequest())
    
    for i, frame in enumerate(frame_iter):
        positions = frame_to_ndarray(frame)

        render_frame(stdscr, positions)

        c = stdscr.getch()

        if c == ord("q"):
            break

    stdscr.clear()
    stdscr.addstr(0, 0, "Closing...")
    stdscr.refresh()

    channel.close()


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
        .frame.arrays[POSITIONS]
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
    arguments = parser.parse_args()
    return arguments


def main(stdscr):
    arguments = handle_user_args()
    write_trajectory_from_server(
        stdscr,
        address=arguments.host,
        port=arguments.port,
    )


if __name__ == '__main__':
    curses.wrapper(main)

import time

screen = curses.initscr()
screen.clear()

time.sleep(3)
