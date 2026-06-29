import numpy as np

def intersects_sphere_line(center: list[float], radius: float, positions: list):
    return np.any(np.sum(np.square(np.subtract(center, positions)), axis=1)) <= radius * radius
