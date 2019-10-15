# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.
"""
Connects to a Narupa simulation and, assuming the simulated molecules to be the
narupa nanotube demo (nanotube with no hydrogens, single methane), computes and
sends over OSC interesting metrics for sonification.
Run with:

.. code bash
    narupa-omm-ase nanotube.xml
    python nanotube.py --osc-port 9000
"""
from osc_app import OscApp

import topology
import itertools
import math

from numpy import arccos, array, dot, pi, cross
from numpy.linalg import norm


# from: https://gist.github.com/nim65s/5e9902cd67f094ce65b0
def distance_segment_point(a, b, p):
    """ segment line AB, point P, where each one is an array([x, y]) """
    if all(a == p) or all(b == p):
        return 0
    if arccos(dot((p - a) / norm(p - a), (b - a) / norm(b - a))) > pi / 2:
        return norm(p - a)
    if arccos(dot((p - b) / norm(p - b), (a - b) / norm(a - b))) > pi / 2:
        return norm(p - b)
    return norm(cross(a - b, a - p)) / norm(b - a)


def clamp01(value):
    return min(1, max(0, value))


def inv_lerp(min, max, value):
    return clamp01((value - min) / (max - min))


def lerp(min, max, u):
    u = clamp01(u)
    return min * (1 - u) + max * u


def centroid(*positions):
    x = y = z = 0
    factor = 1 / len(positions)
    for position in positions:
        x += position[0]
        y += position[1]
        z += position[2]
    return x * factor, y * factor, z * factor


def distance(position_a, position_b):
    dx = position_b[0] - position_a[0]
    dy = position_b[1] - position_a[1]
    dz = position_b[2] - position_a[2]

    return math.sqrt(dx * dx + dy * dy + dz * dz)


def build_frame_generator(osc_client):
    first_frame = osc_client.frame_client.wait_until_first_frame()
    atom_listing = topology.atom_listing_from_frame(first_frame)

    # find the methane's single carbon, knowing it to be the single neighbour of
    # all hydrogens present.
    hydrogen_atoms = [atom for atom in atom_listing if atom['element'] == 1]
    methane_atom = hydrogen_atoms[0]['neighbours'][0]
    methane_index = methane_atom['index']

    # find the nanotube "ends", defined as the terminal carbons with only two
    # neighbours instead of three. Assume two ends lie somewhere on different
    # sides of atom index 30.
    carbons = [atom for atom in atom_listing if atom['element'] == 6]
    ends = [atom for atom in carbons if len(atom['neighbours']) == 2]
    front_indexes = [atom['index'] for atom in ends if atom['index'] < 30]
    back_indexes = [atom['index'] for atom in ends if atom['index'] > 30]

    def frame_to_osc_messages(frame):
        front_positions = [frame.particle_positions[index] for index in front_indexes]
        back_positions = [frame.particle_positions[index] for index in back_indexes]
        front_point = centroid(*front_positions)
        back_point = centroid(*back_positions)

        tube_length = distance(front_point, back_point)
        front_radiuses = [distance(front_point, position) for position in front_positions]
        back_radiuses = [distance(back_point, position) for position in back_positions]

        max_radius = max(itertools.chain(front_radiuses, back_radiuses))
        min_radius = min(itertools.chain(front_radiuses, back_radiuses))

        methane_position = frame.particle_positions[methane_index]
        methane_distance = distance_segment_point(array(front_point),
                                                  array(back_point),
                                                  array(methane_position))

        methane_to_centroid = distance(methane_position,
                                       centroid(front_point, back_point))
        methane_to_centroid_factor = inv_lerp(tube_length / 2, 0, methane_to_centroid)

        interiority = inv_lerp(max_radius, 0, methane_distance)
        centrality = methane_to_centroid_factor * interiority

        yield "/nanotube/length", tube_length
        yield "/nanotube/radius/min", min_radius
        yield "/nanotube/radius/max", max_radius
        yield "/nanotube/radius/ratio", min_radius / max_radius
        yield "/methane/distance", methane_distance
        yield "/methane/interiority", interiority
        yield "/methane/centrality", centrality

        if 'energy.kinetic' in frame.values:
            yield "/energy/kinetic", frame.kinetic_energy

    return frame_to_osc_messages


if __name__ == '__main__':
    app = OscApp(build_frame_generator)
    app.run()
