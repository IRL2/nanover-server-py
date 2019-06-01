"""
Narupa frame client that renders atoms to a curses display in the terminal.

The client connects to a Narupa Frame server, and renders the frames it receives
into the terminal.
"""
import argparse

import math
import numpy as np
import curses

import transformations

import grpc
from narupa.protocol.trajectory import (
    TrajectoryServiceStub,
    GetFrameRequest,
    GetFrameResponse,
)
from narupa.trajectory.frame_data import POSITIONS

cells = [" ", ".", "o", "O", "@"]

def render_positions_to_window(window, positions: np.ndarray):
    xi = 0
    yi = 1
    zi = 2

    h, w = window.getmaxyx()

    positions[:,xi] += (w / 2)
    positions[:,yi] += (h / 2)

    counts = {}
    depths = {}
    indexes = {}

    for index, position in enumerate(positions):
        coord = int(round(position[xi])), int(round(position[yi]))
        count = counts[coord] if coord in counts else 0
        counts[coord] = count + 1

        depth = depths[coord] if coord in depths else -100000

        if position[zi] > depth:
            indexes[coord] = index
        
        depths[coord] = max(depth, position[zi])

    min_depth = min(depths.values())
    max_depth = max(depths.values())

    for coord, count in counts.items():
        x, y = coord

        if x < 0 or x >= w or y < 0 or y >= h:
            continue
        if x == w - 1 and y == h - 1:
            continue

        # by depth
        if True:
            depth = (depths[coord] - min_depth) / (max_depth - min_depth)
            cell_index = int(min(round(depth * len(cells)), len(cells) - 1))
        # by density
        else:
            cell_index = min(math.ceil(count / 3), 4)
        
        index_uniform = (indexes[coord] * 8 / len(positions)) % 1
        color_index = int(min(round(index_uniform * 7), 7 - 1)) + 1

        window.addstr(y, x, cells[cell_index], curses.color_pair(color_index))

def show_controls(window):
    window.addstr(0, 0, "arrow keys -- rotate camera")
    window.addstr(1, 0, " < >  keys -- zoom")

def write_trajectory_from_server(stdscr, *, address: str, port: int):
    """
    Connect to a Narupa frame server and render the received frames to a curses 
    window.

    :param stdscr: curses window to render to.
    :param address: Host name to connect to.
    :param port: Port to connect to on the host.
    """
    stdscr.clear()
    stdscr.nodelay(True)
    stdscr.addstr(0, 0, "Connecting...")
    stdscr.refresh()

    colors = [
        curses.COLOR_BLUE, 
        curses.COLOR_RED,
        curses.COLOR_MAGENTA,
        curses.COLOR_CYAN,
        curses.COLOR_YELLOW,
        curses.COLOR_GREEN,
        curses.COLOR_WHITE
    ]

    for i, color in enumerate(colors):
        curses.init_pair(i + 1, color, curses.COLOR_BLACK)

    host = '{}:{}'.format(address, port)
    channel = grpc.insecure_channel(host)
    stub = TrajectoryServiceStub(channel)
    frame_iter = stub.SubscribeFrames(GetFrameRequest())

    stdscr.clear()
    stdscr.addstr(0, 0, "Connected.")

    angle = 0
    angle2 = 0
    scale = 5

    for i, frame in enumerate(frame_iter):
        positions = frame_to_ndarray(frame)
        
        positions

        center = np.sum(positions, axis=0) / len(positions)
        camera = transformations.scale_matrix(scale, origin=center)
        camera = camera @ transformations.rotation_matrix(angle, [0, 0, 1], point=center)
        camera = camera @ transformations.rotation_matrix(angle2, [1, 0, 0], point=center)

        positions = np.append(positions, np.ones((len(positions), 1)), axis=-1)
        pos = positions[0]
        positions = positions @ camera

        #stdscr.addstr(0, 0, "Frame {}: ({} positions)".format(i, len(positions)))

        stdscr.clear()

        #show_controls(stdscr)

        render_positions_to_window(stdscr, positions)
        show_controls(stdscr)
        stdscr.refresh()

        c = stdscr.getch()

        if c == ord("q"):
            break
        if c == curses.KEY_LEFT:
            angle += 0.1
        if c == curses.KEY_RIGHT:
            angle -= 0.1
        if c == curses.KEY_UP:
            angle2 += 0.1
        if c == curses.KEY_DOWN:
            angle2 -= 0.1
        if c == ord(","):
            scale *= .9
        if c == ord("."):
            scale /= .9

        curses.flushinp()

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
    positions = np.array(raw_positions, dtype=np.float32).reshape((-1, 3))
    return positions


def handle_user_args() -> argparse.Namespace:
    """
    Parse the arguments from the command line.

    :return: The namespace of arguments read from the command line.
    """
    description = "Connect to a Narupa trajectory server and render the atoms live in a curses display."
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
