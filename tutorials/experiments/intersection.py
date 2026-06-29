import numpy as np

def intersects_sphere_line(center: list[float], radius: float, points: list):
    return np.any(np.sum(np.square(np.subtract(center, points)), axis=1)) <= radius * radius
