"""
Narupa frame client that renders atoms to a curses display in the terminal.

The client connects to a Narupa Frame server, and renders the frames it receives
into the terminal.
"""
import argparse

import math
import numpy as np
import curses
import colorsys
import time

import transformations

import grpc
from narupa.protocol.trajectory import (
    TrajectoryServiceStub,
    GetFrameRequest,
    GetFrameResponse,
)
from narupa.trajectory.frame_data import POSITIONS

cells = [" ", ".", "o", "O", "@"]

MAX_COLORS = 16
CHAR_LOOKUP = None
COLOR_LOOKUP = None

def render_positions_to_window(window, positions: np.ndarray):
    xi = 0
    yi = 1
    zi = 2

    h, w = window.getmaxyx()

    positions[:,xi] += (w / 2)
    positions[:,yi] += (h / 2)
    positions[:,zi] *= 1000
    np.round(positions, out=positions)
    positions = positions.astype(int)

    depth_buffer = np.full((w, h), 0, dtype=np.float32)
    index_buffer = np.full((w, h), -1, dtype=np.float32)
    rendered_cells = set()

    min_depth = max_depth = positions[0][zi]

    for index, position in enumerate(positions):
        x, y = position[xi], position[yi]

        if x < 0 or x >= w or y < 0 or y >= h:
            continue
        if x == w - 1 and y == h - 1:
            continue

        coord = y * w + x

        prev_depth = depth_buffer[x, y] if coord in rendered_cells else min_depth
        this_depth = position[zi]

        min_depth = min(min_depth, this_depth)
        max_depth = max(max_depth, this_depth)

        if this_depth >= prev_depth:
            index_buffer[x, y] = index
            depth_buffer[x, y] = this_depth
            rendered_cells.add(coord)

    char_count = len(cells)     
    color_count = MAX_COLORS - 1

    # transform depths into cell indexes
    depth_buffer -= min_depth
    depth_buffer /= (max_depth - min_depth)
    depth_buffer *= len(cells)
    depth_buffer.round(out=depth_buffer)
    depth_buffer = depth_buffer.astype(int)

    # transform particle indexes into color indexes
    index_buffer *= (MAX_COLORS / len(positions))
    index_buffer %= 1
    index_buffer *= color_count
    index_buffer.round(out=index_buffer)
    index_buffer = index_buffer.astype(int)

    char_buffer = CHAR_LOOKUP[depth_buffer]
    color_buffer = COLOR_LOOKUP[index_buffer]

    for y in range(h):
        for x in range(w):
            coord = y * w + x

            if coord not in rendered_cells:
                continue 

            window.addstr(y, x, char_buffer[x, y], color_buffer[x, y])

def show_controls(window):
    window.addstr(0, 0, "arrow keys -- rotate camera")
    window.addstr(1, 0, " < >  keys -- zoom")

def write_trajectory_from_server(stdscr, *, address: str, port: int, custom_colors=False, boxes=False):
    """
    Connect to a Narupa frame server and render the received frames to a curses 
    window.

    :param stdscr: curses window to render to.
    :param address: Host name to connect to.
    :param port: Port to connect to on the host.
    """
    global MAX_COLORS

    stdscr.clear()
    stdscr.nodelay(True)
    stdscr.addstr(0, 0, "Connecting...")
    stdscr.refresh()

    if boxes:
        cells[:] = [" ", "░", "▒", "▓", "█"]

    for i in range(1, MAX_COLORS):
        curses.init_pair(i, i, curses.COLOR_BLACK)

    if curses.can_change_color() and custom_colors:
        curses.init_color(0, 0, 0, 0)
        for i in range(1, MAX_COLORS):
            r, g, b = colorsys.hsv_to_rgb(i / MAX_COLORS, .5, 1)
            curses.init_color(i, int(r * 1000), int(g * 1000), int(b * 1000))
    else:
        MAX_COLORS = 8

    host = '{}:{}'.format(address, port)
    channel = grpc.insecure_channel(host)
    stub = TrajectoryServiceStub(channel)
    frame_iter = stub.SubscribeFrames(GetFrameRequest())

    stdscr.clear()
    stdscr.addstr(0, 0, "Connected.")

    angle = 0
    angle2 = 0
    scale = 5

    global COLOR_LOOKUP, CHAR_LOOKUP
    thing = list(cells)
    thing.append(thing[-1])
    CHAR_LOOKUP = np.array(thing)

    thing = [curses.color_pair(i) for i in range(MAX_COLORS)]
    thing.append(thing[-1])
    COLOR_LOOKUP = np.array(thing)

    for i, frame in enumerate(frame_iter):
        start_time = time.time()

        positions = frame_to_ndarray(frame)

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
        #stdscr.refresh()

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

        h, w = stdscr.getmaxyx()
        stdscr.addstr(h - 1, 0, "{0:.3} fps".format(1.0 / (time.time() - start_time)))
        stdscr.refresh()

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
    parser.add_argument('--rainbow', action="store_true", default=False)
    parser.add_argument('--boxes', action="store_true")
    arguments = parser.parse_args()
    return arguments


def main(stdscr):
    arguments = handle_user_args()
    write_trajectory_from_server(
        stdscr,
        address=arguments.host,
        port=arguments.port,
        custom_colors=arguments.rainbow,
        boxes=arguments.boxes,
    )


if __name__ == '__main__':
    curses.wrapper(main)
