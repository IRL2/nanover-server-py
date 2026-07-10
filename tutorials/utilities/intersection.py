import math
import numpy as np


def closest_point_on_segment(p1, p2, pt):
    p1 = np.array(p1)
    p2 = np.array(p2)
    pt = np.array(pt)

    p2p1 = p2 - p1
    p2p1dot = p2p1.dot(p2p1)

    if p2p1dot == 0:
        return p1, np.linalg.norm(p1 - pt)

    ptp1 = pt - p1

    t = np.dot(ptp1, p2p1) / p2p1dot

    if t < 0:
        closest = p1
    elif t > 1:
        closest = p2
    else:
        closest = p1 + t * p2p1

    return closest, np.linalg.norm(closest - pt)


def closest_point_on_polyline(points: list, target):
    best = (points[0], math.inf, 0)

    for i, (p1, p2) in enumerate(zip(points, points[1:])):
        point, dist = closest_point_on_segment(p1, p2, target)
        if dist < best[1]:
            best = (point, dist, i)

    return best
