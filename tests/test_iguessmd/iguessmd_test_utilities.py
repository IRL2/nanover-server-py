"""
A set of utility objects/functions for testing the iGUESSMD module.
"""

import numpy as np

def define_circular_path(n: int) -> tuple[np.ndarray, np.ndarray]:
    """
    Generate a circular reaction path composed of 2^n + 1 points (where the
    initial and final points overlap) and the associated tangent vectors,
    given the exponent n
    :param n: The exponent n defining the number of points (2^n + 1)
    :return: A tuple (path, tangent_vectors) for the circular path
    """
    # Define no. of points for path
    n_points = 2**n + 1

    # Define unit circular test path (such that final point overlays initial point)
    angles = np.linspace(0, 2 * np.pi, n_points)
    x = np.cos(angles - np.pi / (n_points - 1))
    y = np.sin(angles - np.pi / (n_points - 1))
    z = np.zeros(x.shape)
    circular_path = np.array([x, y, z]).transpose()

    # Define expected unit tangent vectors (including duplicated penultimate tangent)
    x_prime = -np.sin(angles)
    y_prime = np.cos(-angles)
    expected_tangent_vectors = np.array([x_prime, y_prime, z]).transpose()
    expected_tangent_vectors[-1] = expected_tangent_vectors[-2]

    return circular_path, expected_tangent_vectors