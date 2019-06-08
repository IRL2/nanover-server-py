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

from narupa.trajectory import FrameClient

from transformations import rotation_matrix, scale_matrix

character_sets = {
    "boxes": ["░", "▒", "▓", "█"],
    "blobs": [".", "-", "+", "o", "O", "@"],
    "extended-blobs": ["·", "-", "+", "◌", "○", "ø", "●", "■"],
}

character_sets_indexed = list(character_sets.values())

class UserQuitException(Exception):
    pass

class Camera:
    def __init__(self):
        self.angle = 0
        self.pitch = math.pi / 2
        self.scale = 5
        self.origin = None

    def get_matrix(self):
        return (scale_matrix(self.scale, origin=self.origin)
              @ rotation_matrix(self.angle, [0, 0, 1], point=self.origin)
              @ rotation_matrix(self.pitch, [1, 0, 0], point=self.origin))

class Timer:
    def __init__(self):
        self.last_time = time.time()

    def get_delta(self):
        return time.time() - self.last_time

    def reset(self):
        delta_time = self.get_delta()
        self.last_time = time.time()
        return delta_time

class Style:
    def __init__(self, characters, colors):
        self.set_characters(characters)
        self.set_colors(colors)

    def set_characters(self, characters):
        self.characters = characters
        self.character_lookup = np.array(self.characters)

    def set_colors(self, colors):
        self.colors = colors
        self.color_lookup = np.array(self.colors)

class Renderer:
    def __init__(self, window, style):
        self.window = window
        self.style = style
        self.camera = Camera()

        self.frame = None
        self.elements = None

        self.cpk = True

    def render(self):
        self.window.clear()
        
        if self.frame:
            positions = self._process_frame(self.frame)
            render_positions_to_window(self.window, positions, self.style, self.elements if self.cpk else None)
        
        self.window.noutrefresh()

    def _process_frame(self, frame):
        # convert positions from frame
        positions = np.array(self.frame.particle_positions, dtype=np.float32).reshape((-1, 3))

        # center the camera
        self.camera.origin = np.average(positions, axis=0)

        # add w component to positions then multiply by camera matrix
        positions = np.append(positions, np.ones((len(positions), 1)), axis=-1)
        positions = positions @ self.camera.get_matrix()

        return positions

element_colors = {
    #1: curses.COLOR_WHITE,
    6: curses.COLOR_CYAN,
    7: curses.COLOR_BLUE,
    8: curses.COLOR_RED,
    16: curses.COLOR_YELLOW,
}

def render_positions_to_window(window, positions: np.ndarray, renderer, elements=None):
    h, w = window.getmaxyx()

    positions[:,0] += (w / 2)
    positions[:,1] += (h / 2)
    positions[:,2] *= 1000
    np.round(positions, out=positions)
    positions = positions.astype(int)

    depth_buffer = np.full((w, h), 0, dtype=np.float32)
    color_buffer = {}

    color_count = len(renderer.colors) - 1

    min_depth = max_depth = positions[0][2]

    def record_color(index, x, y, z):
        nonlocal min_depth, max_depth

        coord = (x, y)

        prev_depth = depth_buffer[x, y] if coord in color_buffer else min_depth
        this_depth = z

        min_depth = min(min_depth, this_depth)
        max_depth = max(max_depth, this_depth)

        if this_depth >= prev_depth:
            if elements:
                element = elements[index]
                if element not in element_colors:
                    return

                color_index = element_colors[element]
            else:
                color_index = int((index * .1) % color_count)

            color_buffer[coord] = renderer.colors[color_index]
            depth_buffer[x, y] = this_depth

    for index, position in enumerate(positions):
        x, y, z = position[0], position[1], position[2]

        if x < 0 or x >= w or y < 0 or y >= h:
            continue
        if x == w - 1 and y == h - 1:
            continue

        record_color(index, x, y, z)

    char_count = len(renderer.characters)
    depth_scale = char_count / (max_depth - min_depth)

    # transform depths into cell indexes
    depth_buffer -= min_depth
    depth_buffer *= depth_scale
    np.rint(depth_buffer, out=depth_buffer)

    depth_buffer.clip(0, char_count-1, out=depth_buffer)
    depth_buffer = depth_buffer.astype(int)

    char_buffer = renderer.character_lookup[depth_buffer]

    for (x, y), color in color_buffer.items():
        window.addch(y, x, char_buffer[x, y], color)

def generate_colors(override_colors=False):
    max_colors = 8

    if curses.can_change_color() and override_colors:
        max_colors = 16
        curses.init_color(0, 0, 0, 0)
        for i in range(1, max_colors):
            r, g, b = colorsys.hsv_to_rgb(i / max_colors, .5, 1)
            curses.init_color(i, int(r * 1000), int(g * 1000), int(b * 1000))

    for i in range(1, max_colors):
        curses.init_pair(i, i, curses.COLOR_BLACK)

    colors = [curses.color_pair(i) for i in range(max_colors)]

    return colors

def run_curses_client(stdscr, *, address: str, port: int, override_colors=False, skin="boxes", cpk=False):
    """
    Connect to a Narupa frame server and render the received frames to a curses 
    window.

    :param stdscr: Curses window to render to.
    :param address: Host name to connect to.
    :param port: Port to connect to on the host.
    :param override_colors: Whether to modify the terminal colors.
    :param skin: Name of character set to use.
    """

    colors = generate_colors(override_colors=override_colors)
    style = Style(character_sets[skin], colors)

    fps_timer = Timer()
    render_timer = Timer()
    renderer = Renderer(stdscr, style)

    curses.curs_set(False)
    stdscr.clear()
    stdscr.nodelay(True)
    stdscr.addstr(0, 0, "Connecting...")
    stdscr.refresh()

    client = FrameClient(address=address, port=port)

    stdscr.clear()
    stdscr.addstr(0, 0, "Connected.")

    def get_frame(frame_index, frame):
        renderer.frame = frame

        if 'particle.element' in frame:
            renderer.elements = frame.particle_elements

    client.subscribe_frames_async(get_frame)

    def show_controls(window):
        window.addstr(0, 0, "arrows -- rotate camera")
        window.addstr(1, 0, " < >   -- zoom")
        window.addstr(2, 0, "  x    -- toggle cpk colors")
        window.addstr(3, 0, "  z    -- cycle skins")

    def rotate_plus():
        renderer.camera.angle += .1
    def rotate_minus():
        renderer.camera.angle -= .1
    def pitch_plus():
        renderer.camera.pitch += .1
    def pitch_minus():
        renderer.camera.pitch -= .1
    def zoom_in():
        renderer.camera.scale *= .9
    def zoom_out():
        renderer.camera.scale /= .9
    def toggle_cpk():
        renderer.cpk = not renderer.cpk
    def cycle_charsets():
        index = character_sets_indexed.index(renderer.style.characters)
        index = (index + 1) % len(character_sets_indexed)
        renderer.style.set_characters(character_sets_indexed[index])
    def quit():
        raise UserQuitException()

    bindings = {
        curses.KEY_LEFT:  rotate_plus,
        curses.KEY_RIGHT: rotate_minus,
        curses.KEY_UP:    pitch_plus,
        curses.KEY_DOWN:  pitch_minus,
        ord(","):         zoom_in,
        ord("."):         zoom_out,
        ord("q"):         quit,
        ord("x"):         toggle_cpk,
        ord("z"):         cycle_charsets,
    }

    def check_input():
        char = stdscr.getch()

        if char in bindings:
            bindings[char]()

        curses.flushinp()

    h, w = stdscr.getmaxyx()

    def loop():
        check_input()

        if render_timer.get_delta() < 0.05:
            return

        render_timer.reset()

        stdscr.clear()

        renderer.render()
        show_controls(stdscr)
        stdscr.addstr(h - 1, 0, "{0:.3} fps".format(1.0 / fps_timer.reset()))
        
        stdscr.noutrefresh()
        curses.doupdate()

    try:
        while True:
            loop()
            time.sleep(.01)
    except UserQuitException:
        pass
    finally:
        stdscr.clear()
        stdscr.addstr(0, 0, "Closing...")
        stdscr.refresh()
        client.close()


def handle_user_args() -> argparse.Namespace:
    """
    Parse the arguments from the command line.

    :return: The namespace of arguments read from the command line.
    """
    description = "Connect to a Narupa trajectory server and render the atoms live in a curses display."
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('--host', default='localhost')
    parser.add_argument('--port', type=int, default=8000)
    parser.add_argument('--rainbow', action="store_true")
    parser.add_argument('--skin', default="boxes")
    parser.add_argument('--cpk', action="store_true")
    arguments = parser.parse_args()
    return arguments


def main(stdscr):
    arguments = handle_user_args()
    run_curses_client(
        stdscr,
        address=arguments.host,
        port=arguments.port,
        override_colors=arguments.rainbow,
        skin=arguments.skin,
        cpk=arguments.cpk,
    )


if __name__ == '__main__':
    curses.wrapper(main)
