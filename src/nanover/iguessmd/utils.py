import numpy as np


def get_every_nth(
    array: np.ndarray, axis: int, every_nth: int, include_end: bool = True
) -> np.ndarray:
    """
    Returns a reduced version of the input array composed of every nth
    value of the original array along a defined axis. By default the
    final entry of the axis is also included, regardless of the value
    of the stride (every_nth).
    :param array: Original array to be reduced
    :param axis: Axis along which to reduce the array
    :param every_nth: (int | None) Defines the stride for which to retrieve
      every nth point of the array
    :param include_end: Bool defining whether to include the last entry of the array
    :return: Reduced version of original array
    """
    # Assert array has the dimensions required
    assert len(array.shape) - 1 >= axis

    # Throw error if requested stride is zero
    if every_nth <= 0:
        raise ValueError("Every nth value must be a non-zero positive integer")

    # Return original array if requested stride is 1
    if every_nth == 1:
        return array

    # Define which entries of the array to take
    indices = np.arange(0, array.shape[axis], every_nth)

    # Optionally include the final element if not automatically included
    if (array.shape[axis] - 1) % every_nth != 0 and include_end:
        indices = np.append(indices, array.shape[axis] - 1)

    return np.take(array, indices, axis=axis)


def calculate_unit_tangents(
    path: np.ndarray,
) -> np.ndarray:
    """
    Calculate the unit tangents of a path in 3D space using the
    forward difference approximation, setting the final tangent
    equal to the penultimate tangent.
    :param path: (k * 3) array of k points defining the path
    :return: (k * 3) array of unit tangents to the path
    """
    # Initialise tangents array
    unit_tangents = np.zeros(path.shape)

    # Calculate displacements between path points
    displacements = np.diff(path, axis=0)

    # Normalise displacement vectors
    displacements /= np.linalg.norm(displacements, axis=1, keepdims=True)

    # Populate unit tangents
    unit_tangents[:-1] = displacements
    unit_tangents[-1] = displacements[-1]

    return unit_tangents
