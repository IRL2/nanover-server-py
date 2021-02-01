# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.
"""
Narupa frame client that renders atoms to a curses display in the terminal.

The client connects to a Narupa Frame server, and renders the frames it receives
into the terminal.
"""
import sys
import textwrap
from typing import Dict, Callable, Sequence, Any

from narupa.app.app_server import DEFAULT_NARUPA_PORT

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

from narupa.app import NarupaImdClient

from transformations import rotation_matrix, scale_matrix
import rendering


class Camera:
    def __init__(self):
        self.angle = 0
        self.pitch = math.pi / 2
        self.scale = 5

    def get_matrix(self) -> np.ndarray:
        return (scale_matrix(self.scale)
                @ rotation_matrix(self.angle, [0, 0, 1])
                @ rotation_matrix(self.pitch, [1, 0, 0]))


class FPSTimer:
    """
    Stores a number of frame intervals to compute an average frames per second
    measure.
    """
    intervals: Sequence[float]
    fps: float
    count: int

    def __init__(self, count: int = 5):
        self.intervals = []
        self.fps = 0
        self.count = count
        self._prev_checkpoint = time.monotonic()

    def checkpoint(self):
        interval = time.monotonic() - self._prev_checkpoint
        self._prev_checkpoint = time.monotonic()
        self.intervals.append(interval)
        self.intervals[:] = self.intervals[-self.count:]
        try:
            self.fps = len(self.intervals) / sum(self.intervals)
        except ZeroDivisionError:
            pass


class FrameTimer:
    def __init__(self, fps: int = 30):
        self.interval = 1 / fps
        self.next_time = time.monotonic()

    def check(self) -> bool:
        """
        Return True if the interval has passed, resetting the timer for the next
        interval if it has.
        """
        now = time.monotonic()
        if now >= self.next_time:
            self.next_time = now + self.interval
            return True
        return False


class Renderer:
    charset: Sequence[str]
    colors: Sequence[rendering.ColorPair]
    shader: rendering.Shader
    camera: Camera

    def __init__(self,
                 window,
                 charset: Sequence[str],
                 colors: Sequence[rendering.ColorPair],
                 shader: rendering.Shader):
        self.window = window
        self.charset = charset
        self.colors = colors
        self.shader = shader
        self.camera = Camera()

        self.positions = None
        self.elements = None
        self.bonds = None

        self.show_hydrogens = True

        self._offset = (0, 0, 0)
        self._recenter = True

    def render(self):
        self.window.clear()
        
        if self.positions is not None:
            self._process_frame()

            frame = {
                'positions': self.positions,
                'elements': self.elements,
                'bonds': self.bonds,
                'skip_atoms': set(),
            }

            if not self.show_hydrogens:
                frame['skip_atoms'] = set(index for index, element in enumerate(frame['elements']) if element == 1)

            pixels = self.shader(frame, len(self.colors) - 1)
            rendering.render_pixels_to_window(self.window, self.charset,
                                              self.colors, pixels)
        
        self.window.noutrefresh()

    def recenter_camera(self):
        self._recenter = True

    def _process_frame(self):
        # center the camera
        if self._recenter:
            self._offset = np.median(self.positions, axis=0)
            self._recenter = False
            
        self.positions -= self._offset

        # add w component to positions then multiply by camera matrix
        self.positions = np.append(self.positions, np.ones((len(self.positions), 1)), axis=-1)
        self.positions = self.positions @ self.camera.get_matrix()

        # move coordinate origin to window center
        h, w = self.window.getmaxyx()
        self.positions[:, 0] += (w / 2)
        self.positions[:, 1] += (h / 2)


class CursesFrontend:
    client: NarupaImdClient
    renderer: Renderer
    bindings: Dict[int, Callable]

    def __init__(self, stdscr, client: NarupaImdClient, override_colors=False):
        self.stdscr = stdscr
        self.client = client
        self.user_quit = False

        self._render_timer = FrameTimer(fps=30)
        self._fps_timer = FPSTimer(5)

        colors = generate_colors(override_colors=override_colors)
        self.renderer = Renderer(stdscr,
                                 rendering.character_sets["boxes"],
                                 colors,
                                 rendering.SHADERS[0])

        self.bindings = {
            curses.KEY_LEFT: self.rotate_camera_clockwise,
            curses.KEY_RIGHT: self.rotate_camera_anticlockwise,
            curses.KEY_UP: self.pitch_camera_up,
            curses.KEY_DOWN: self.pitch_camera_down,
            ord(","): self.zoom_camera_in,
            ord("."): self.zoom_camera_out,
            ord("q"): self.quit,
            ord("x"): self.cycle_shaders,
            ord("z"): self.cycle_charsets,
            ord("c"): self.renderer.recenter_camera,
            ord("v"): self.toggle_hydrogens,
        }

    def run(self):
        curses.curs_set(False)
        self.stdscr.clear()
        self.stdscr.nodelay(True)
        self.stdscr.addstr(0, 0, "Starting...")
        self.stdscr.refresh()

        while not self.user_quit:
            self.update()
            time.sleep(.001)

    def update(self):
        self.check_input()

        if not self.client.first_frame:
            return

        if self._render_timer.check():
            self.render()

    def render(self):
        self.stdscr.clear()
        self.renderer.positions = np.array(self.client.latest_frame.particle_positions, dtype=np.float32)
        self.renderer.bonds = self.client.first_frame.bond_pairs
        self.renderer.elements = self.client.first_frame.particle_elements

        self.renderer.render()
        self.show_controls()

        self._fps_timer.checkpoint()

        h, w = self.stdscr.getmaxyx()
        self.stdscr.addstr(h - 1, 0, "{} fps".format(round(self._fps_timer.fps)))
        self.stdscr.noutrefresh()
        curses.doupdate()

    def show_controls(self):
        self.stdscr.addstr(0, 0, "arrows -- rotate camera")
        self.stdscr.addstr(1, 0, "  c    -- recenter camera")
        self.stdscr.addstr(2, 0, " < >   -- zoom")

        self.stdscr.addstr(4, 0, " x -- cycle representations")
        self.stdscr.addstr(5, 0, " z -- cycle skins")
        self.stdscr.addstr(6, 0, " v -- toggle hydrogens")
        self.stdscr.addstr(7, 0, " q -- quit")

    def check_input(self):
        char = self.stdscr.getch()

        if char in self.bindings:
            self.bindings[char]()

        curses.flushinp()

    def rotate_camera_clockwise(self):
        self.renderer.camera.angle += .1

    def rotate_camera_anticlockwise(self):
        self.renderer.camera.angle -= .1

    def pitch_camera_up(self):
        self.renderer.camera.pitch += .1

    def pitch_camera_down(self):
        self.renderer.camera.pitch -= .1

    def zoom_camera_in(self):
        self.renderer.camera.scale *= .9

    def zoom_camera_out(self):
        self.renderer.camera.scale /= .9

    def cycle_shaders(self):
        self.renderer.shader = next_in_cycle(rendering.SHADERS,
                                             self.renderer.shader)

    def cycle_charsets(self):
        charset = next_in_cycle(rendering.character_sets_indexed,
                                self.renderer.charset)
        self.renderer.charset = charset

    def toggle_hydrogens(self):
        self.renderer.show_hydrogens = not self.renderer.show_hydrogens

    def quit(self):
        self.user_quit = True


def override_terminal_colors(max_colors: int,
                             saturation: float = .5,
                             value: float = 1):
    """
    Updates the terminal colors to black + a hue cycle.
    """
    curses.init_color(0, 0, 0, 0)
    for i in range(1, max_colors):
        r, g, b = colorsys.hsv_to_rgb(i / max_colors, saturation, value)
        curses.init_color(i, int(r * 1000), int(g * 1000), int(b * 1000))


def generate_colors(override_colors=False) -> Sequence[rendering.ColorPair]:
    """
    Return a sequence of all available curses color pairs of a foreground color
    on a black background.
    """
    max_colors = 8

    if curses.can_change_color() and override_colors:
        max_colors = 16
        override_terminal_colors(max_colors)

    for i in range(1, max_colors):
        curses.init_pair(i, i, curses.COLOR_BLACK)

    colors = [curses.color_pair(i) for i in range(max_colors)]

    return colors


def next_in_cycle(cycle: Sequence[Any], current: Any, move: int = 1):
    current_index = cycle.index(current)
    next_index = (current_index + move) % len(cycle)
    return cycle[next_index]


def handle_user_args(args=None) -> argparse.Namespace:
    """
    Parse the arguments from the command line.

    :return: The namespace of arguments read from the command line.
    """
    description = textwrap.dedent("""\
    Connect to a Narupa trajectory server and render the atoms live in a curses 
    display.
    """)
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('--autoconnect', nargs='?', default=False, metavar='NAME')
    parser.add_argument('--hostname', default=None)
    parser.add_argument('--port', '-p', default=None, type=int)
    parser.add_argument('--rainbow', '-r', action="store_true")
    arguments = parser.parse_args(args)
    return arguments


def main(stdscr):
    arguments = handle_user_args()

    if arguments.autoconnect is not False:
        client = NarupaImdClient.autoconnect(name=arguments.autoconnect)
    else:
        address = (
            arguments.hostname or 'localhost',
            arguments.port or DEFAULT_NARUPA_PORT,
        )
        client = NarupaImdClient(trajectory_address=address)

    with client:
        telmol = CursesFrontend(stdscr, client, override_colors=arguments.rainbow)
        telmol.run()


if __name__ == '__main__':
    curses.wrapper(main)
