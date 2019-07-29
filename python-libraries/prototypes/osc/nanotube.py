from osc_app import OscApp

import selections
import itertools
import math

from numpy import arccos, array, dot, pi, cross
from numpy.linalg import norm


# from: https://gist.github.com/nim65s/5e9902cd67f094ce65b0
def distance_segment_point(A, B, P):
    """ segment line AB, point P, where each one is an array([x, y]) """
    if all(A == P) or all(B == P):
        return 0
    if arccos(dot((P - A) / norm(P - A), (B - A) / norm(B - A))) > pi / 2:
        return norm(P - A)
    if arccos(dot((P - B) / norm(P - B), (A - B) / norm(A - B))) > pi / 2:
        return norm(P - B)
    return norm(cross(A-B, A-P))/norm(B-A)


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
    atom_listing = selections.atom_listing_from_frame(first_frame)

    hydrogen_atoms = [atom for atom in atom_listing if atom['element'] == 1]
    carbon_atom = hydrogen_atoms[0]['neighbours'][0]

    carbon_index = carbon_atom['index']

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

        methane_position = frame.particle_positions[carbon_index]
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
        yield "/energy/kinetic", frame.values['energy.kinetic']

    osc_client.message_generator = frame_to_osc_messages


if __name__ == '__main__':
    app = OscApp(build_frame_generator)
    app.run()
