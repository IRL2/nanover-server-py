"""
Narupa frame client that renders atoms to a curses display in the terminal.

The client connects to a Narupa Frame server, and renders the frames it receives
into the terminal.
"""
import sys

try:
    import curses
except ModuleNotFoundError:
    print("Unable to import the curses module. Try `pip install windows-curses` if you are on windows.")
    sys.exit(1)

import argparse

import math
import numpy as np
import colorsys
import time
import itertools

from narupa.core import DEFAULT_CONNECT_ADDRESS
from narupa.trajectory.frame_server import DEFAULT_PORT
from narupa.trajectory import FrameClient
from narupa.trajectory.frame_data import POSITIONS, ELEMENTS, BONDS

from transformations import rotation_matrix, scale_matrix, translation_matrix, compose_matrix
from bresenham import get_line

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
    def __init__(self, window, style, shader):
        self.window = window
        self.style = style
        self.shader = shader
        self.camera = Camera()

        self.frame = None
        self.elements = None

        self.show_hydrogens = True

        self._offset = (0, 0, 0)
        self._recenter = False

    def render(self):
        self.window.clear()
        
        if self.frame:
            positions = self._process_frame(self.frame)

            frame = dict(positions=positions,
                         elements=self.elements,
                         bonds=self.bonds,
                         skip_atoms=set())

            if not self.show_hydrogens:
                frame['skip_atoms'] = set(index for index, element in enumerate(frame['elements']) if element == 1)

            render_frame_to_window(self.window, self.style, frame, self.shader)
        
        self.window.noutrefresh()

    def recenter_camera(self):
        self._recenter = True

    def _process_frame(self, frame):
        # convert positions from frame
        raw_positions = frame._raw.arrays[POSITIONS].float_values.values
        positions = np.array(raw_positions, dtype=np.float32).reshape((-1, 3))

        # center the camera
        #self.camera.origin = np.median(positions, axis=0)
                
        if self._recenter:
            self._offset = np.median(positions, axis=0)[:3]
            self._recenter = False
            
        positions -= self._offset

        # add w component to positions then multiply by camera matrix
        positions = np.append(positions, np.ones((len(positions), 1)), axis=-1)
        positions = positions @ self.camera.get_matrix()

        center_positions_in_window(self.window, positions)

        return positions

element_colors = {
    1: curses.COLOR_WHITE,
    6: curses.COLOR_CYAN,
    7: curses.COLOR_BLUE,
    8: curses.COLOR_RED,
    16: curses.COLOR_YELLOW,
}

depth_buffer = None

def iterate_frame_atoms(frame):
    for index, position, element in zip(itertools.count(), frame['positions'], frame['elements']):
        if index not in frame['skip_atoms']:
            yield index, position, element

def iterate_frame_bonds(frame):
    for atom_a_index, atom_b_index in frame['bonds']:
        if atom_a_index not in frame['skip_atoms'] and atom_b_index not in frame['skip_atoms']:
            yield atom_a_index, atom_b_index

def atoms_to_pixels_cpk(frame, color_count):
    for index, position, element in iterate_frame_atoms(frame):        
        if element not in element_colors:
            continue

        x, y, z = int(position[0]), int(position[1]), position[2]

        yield element_colors[element], x, y, z

def atoms_to_pixels_gradient(frame, color_count):
    for index, position, element in iterate_frame_atoms(frame):
        x, y, z = int(position[0]), int(position[1]), position[2]
        color_index = int((index * .1) % color_count)

        yield color_index, x, y, z

def bonds_to_pixels_gradient(frame, color_count):
    for atom_a_index, atom_b_index in iterate_frame_bonds(frame):
        start = frame['positions'][atom_a_index]
        end = frame['positions'][atom_b_index]

        start = (int(start[0]), int(start[1]), start[2])
        end = (int(end[0]), int(end[1]), end[2])

        color_index = int(((atom_a_index + atom_b_index) *.05) % color_count)

        for x, y, z in get_line(start, end):
            yield color_index, x, y, z

def render_pixels_to_window(window, style, pixels):
    global depth_buffer
    h, w = window.getmaxyx()

    if depth_buffer is None: 
        depth_buffer = np.full((w, h), 0, dtype=np.float32)
    else:
        depth_buffer.fill(0)

    color_buffer = {}

    minus_infinity = float("-inf")

    def write_pixel(color, x, y, z):
        coord = (x, y)

        prev_depth = depth_buffer[x, y] if coord in color_buffer else minus_infinity
        this_depth = z

        if this_depth >= prev_depth:
            color_buffer[coord] = color
            depth_buffer[x, y] = this_depth

    for color_index, x, y, z in pixels:
        if x < 0 or x >= w or y < 0 or y >= h:
            continue
        if x == w - 1 and y == h - 1:
            continue

        write_pixel(style.colors[color_index], x, y, z)

    min_depth = depth_buffer.min()
    max_depth = depth_buffer.max()

    if not color_buffer:
        return

    char_count = len(style.characters)
    depth_scale = char_count / (max_depth - min_depth)

    # transform depths into cell indexes
    depth_buffer -= min_depth
    depth_buffer *= depth_scale
    np.rint(depth_buffer, out=depth_buffer)
    depth_buffer.clip(0, char_count-1, out=depth_buffer)

    char_buffer = style.character_lookup[depth_buffer.astype(int)]

    for (x, y), color in color_buffer.items():
        window.addch(y, x, char_buffer[x, y], color)

def center_positions_in_window(window, positions):
    h, w = window.getmaxyx()
    positions[:,0] += (w / 2)
    positions[:,1] += (h / 2)

def render_frame_to_window(window, style, frame, shader):
    pixels = shader(frame, len(style.colors) - 1)
    render_pixels_to_window(window, style, pixels)

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

def run_curses_client(stdscr, *, address: str, port: int, override_colors=False):
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
    style = Style(character_sets["boxes"], colors)    
    shaders = [atoms_to_pixels_cpk, atoms_to_pixels_gradient, bonds_to_pixels_gradient]

    fps_timer = Timer()
    render_timer = Timer()
    renderer = Renderer(stdscr, style, shaders[0])

    fpses = []

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

        if BONDS in frame:
            renderer.bonds = frame.bonds
        if ELEMENTS in frame:
            renderer.elements = frame.particle_elements

    client.subscribe_last_frames_async(get_frame)

    def show_controls(window):
        window.addstr(0, 0, "arrows -- rotate camera")
        window.addstr(1, 0, "  c    -- recenter camera")
        window.addstr(2, 0, " < >   -- zoom")

        window.addstr(4, 0, " x -- cycle representations")
        window.addstr(5, 0, " z -- cycle skins")
        window.addstr(6, 0, " v -- toggle hydrogens")

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
    def cycle_shaders():
        index = shaders.index(renderer.shader)
        index = (index + 1) % len(shaders)
        renderer.shader = shaders[index]
    def cycle_charsets():
        index = character_sets_indexed.index(renderer.style.characters)
        index = (index + 1) % len(character_sets_indexed)
        renderer.style.set_characters(character_sets_indexed[index])
    def toggle_hydrogens():
        renderer.show_hydrogens = not renderer.show_hydrogens
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
        ord("x"):         cycle_shaders,
        ord("z"):         cycle_charsets,
        ord("c"):         renderer.recenter_camera,
        ord("v"):         toggle_hydrogens,
    }

    def check_input():
        char = stdscr.getch()

        if char in bindings:
            bindings[char]()

        curses.flushinp()

    h, w = stdscr.getmaxyx()

    frame_delay = 1 / 30

    def loop():
        check_input()

        if render_timer.get_delta() < frame_delay:
            return

        render_timer.reset()

        stdscr.clear()

        renderer.render()
        show_controls(stdscr)

        fpses.append(fps_timer.reset())
        
        if len(fpses) > 5:
            del fpses[0]

        average_fps = sum(fpses) / len(fpses)

        try:
            stdscr.addstr(h - 1, 0, "{} fps".format(round(1.0 / average_fps)))
        except ZeroDivisionError:
            pass

        stdscr.noutrefresh()
        curses.doupdate()

    def main_loop():
        try:
            while True:
                loop()
                time.sleep(.0001)
        except UserQuitException:
            pass
        finally:
            client.close()
    
    main_loop()

def handle_user_args() -> argparse.Namespace:
    """
    Parse the arguments from the command line.

    :return: The namespace of arguments read from the command line.
    """
    description = "Connect to a Narupa trajectory server and render the atoms live in a curses display."
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('--host', default=DEFAULT_CONNECT_ADDRESS)
    parser.add_argument('--port', type=int, default=DEFAULT_PORT)
    parser.add_argument('--rainbow', action="store_true")
    arguments = parser.parse_args()
    return arguments


def main(stdscr):
    arguments = handle_user_args()
    run_curses_client(stdscr,
                      address=arguments.host,
                      port=arguments.port,
                      override_colors=arguments.rainbow)

if __name__ == '__main__':
    curses.wrapper(main)
