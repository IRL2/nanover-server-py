import numpy as np


def intersects_sphere_line(center: list[float], radius: float, positions: list):
    length2 = np.sum(np.subtract(center, positions) ** 2, axis=1)
    return np.any(length2 <= radius * radius)
