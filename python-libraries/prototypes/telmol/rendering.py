# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.
import curses
import itertools
from typing import Tuple, Iterable, Iterator, Callable, Sequence, Dict, Any, Protocol
import numpy as np
from bresenham import get_line

character_sets = {
    "boxes": ["░", "▒", "▓", "█"],
    "blobs": [".", "-", "+", "o", "O", "@"],
    "extended-blobs": ["·", "-", "+", "◌", "○", "ø", "●", "■"],
}

character_sets_indexed = list(character_sets.values())


element_colors = {
    1: curses.COLOR_WHITE,
    6: curses.COLOR_CYAN,
    7: curses.COLOR_BLUE,
    8: curses.COLOR_RED,
    16: curses.COLOR_YELLOW,
}

ColorPair = Any # curses.color_pair
Position = Tuple[float, float, float, float]
Fragment = Tuple[int, int, int, float]
BondPair = Tuple[int, int]
AtomInfo = Tuple[int, Position, int]


# Mypy does not allow to assign to a method. Though, it is confused and thinks
# any attribute declared as a callable is a method. This cause issues in the
# Renreder class in curse_client.py. This is a workaround described in
# https://github.com/python/mypy/issues/708#issuecomment-667989040
class Shader(Protocol):
    def __call__(self, frame: Dict, color_count: int) -> Iterator[Fragment]:
        ...


def iterate_frame_atoms(frame) -> Iterator[AtomInfo]:
    for index, position, element in zip(itertools.count(), frame['positions'], frame['elements']):
        if index not in frame['skip_atoms']:
            yield index, position, element


def iterate_frame_bonds(frame) -> Iterator[BondPair]:
    for atom_a_index, atom_b_index in frame['bonds']:
        if atom_a_index not in frame['skip_atoms'] and atom_b_index not in frame['skip_atoms']:
            yield atom_a_index, atom_b_index


def atoms_to_pixels_cpk(frame, color_count) -> Iterator[Fragment]:
    for index, position, element in iterate_frame_atoms(frame):
        if element not in element_colors:
            continue

        x, y, z = int(position[0]), int(position[1]), position[2]

        yield element_colors[element], x, y, z


def atoms_to_pixels_gradient(frame, color_count) -> Iterator[Fragment]:
    for index, position, element in iterate_frame_atoms(frame):
        x, y, z = int(position[0]), int(position[1]), position[2]
        color_index = int((index * .1) % color_count)

        yield color_index, x, y, z


def bonds_to_pixels_gradient(frame, color_count) -> Iterator[Fragment]:
    for atom_a_index, atom_b_index in iterate_frame_bonds(frame):
        start = frame['positions'][atom_a_index]
        end = frame['positions'][atom_b_index]

        start = (int(start[0]), int(start[1]), start[2])
        end = (int(end[0]), int(end[1]), end[2])

        color_index = int(((atom_a_index + atom_b_index) *.05) % color_count)

        for x, y, z in get_line(start, end):
            yield color_index, x, y, z


def render_pixels_to_window(window,
                            charset: Sequence[str],
                            colors: Sequence[ColorPair],
                            pixels: Iterable[Fragment]):
    h, w = window.getmaxyx()
    depth_buffer = np.full((w, h), 0, dtype=np.float32)
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

        write_pixel(colors[color_index], x, y, z)

    min_depth = depth_buffer.min()
    max_depth = depth_buffer.max()

    if not color_buffer:
        return

    char_count = len(charset)
    depth_scale = char_count / (max_depth - min_depth)

    # transform depths into cell indexes
    depth_buffer -= min_depth
    depth_buffer *= depth_scale
    np.rint(depth_buffer, out=depth_buffer)
    depth_buffer.clip(0, char_count-1, out=depth_buffer)

    character_lookup = np.array(charset)
    char_buffer = character_lookup[depth_buffer.astype(int)]

    for (x, y), color in color_buffer.items():
        window.addch(y, x, char_buffer[x, y], color)


SHADERS = [atoms_to_pixels_cpk,
           atoms_to_pixels_gradient,
           bonds_to_pixels_gradient]
