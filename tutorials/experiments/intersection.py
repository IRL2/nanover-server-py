import numpy as np


def intersects_sphere_line(center: list[float], radius: float, positions: list):
    deltas = np.subtract(center, positions)
    length2 = np.sum(np.square(deltas), axis=1)
    return np.any(length2 <= radius * radius)
