# Utilising some code from https://github.com/marrink-lab/vermouth-martinize.
#
# Copyright 2018 University of Groningen
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""
Tests minimum image convention and distance calculations.
"""

import numpy as np
import pytest
from hypothesis import strategies, given, example

from narupa.imd.imd_force import _calculate_diff_and_sqr_distance


@strategies.composite
def vector_with_random_distance(draw):
    """
    Generate a vector with a random length and orientation.
    The vector is returned as the two points at its extremity.
    Returns
    -------
    point1: numpy.ndarray
    point2: numpy.ndarray
    distance: float
    """
    # Generate random polar coordinates and convert them to euclidean
    # coordinates.
    length = strategies.floats(min_value=0, max_value=100,
                               allow_nan=False, allow_infinity=False)
    angle = strategies.floats(min_value=0, max_value=2 * np.pi,
                              allow_nan=False, allow_infinity=False)
    distance = draw(length)
    theta = draw(angle)
    phi = draw(angle)
    shift = np.array([draw(length), draw(length), draw(length)])
    x = distance * np.sin(theta) * np.cos(phi)  # pylint: disable=invalid-name
    y = distance * np.sin(theta) * np.sin(phi)  # pylint: disable=invalid-name
    z = distance * np.cos(theta)  # pylint: disable=invalid-name
    return shift, np.array([x, y, z]) + shift, distance


@given(vector_with_random_distance())
@example((np.zeros((3,)), np.zeros((3,)), 0))
def test_difference_non_periodic(vector_distance):
    """
    Tests difference and squared distance calculation in non periodic case.
    """
    point1, point2, distance = vector_distance
    diff, dist_sqr = _calculate_diff_and_sqr_distance(point1, point2)
    assert dist_sqr == pytest.approx(distance * distance)
    expected_diff = point1 - point2
    assert np.allclose(diff, expected_diff)


@strategies.composite
def points_in_periodic_box(draw):
    """
    Generates a random periodic box, and two random positions.
    The two points, the periodic box, and the difference and square distance between the two points under
    minimum image convention are returned.
    """

    # Generate random lengths to get a periodic box.
    # box length has to be nonzero.
    length = strategies.floats(min_value=1e-12, max_value=100,
                               allow_nan=False, allow_infinity=False)

    periodic_box_lengths = np.array([draw(length) for x in range(3)])

    # pick two random points in lowest quadrant of the box.
    lengths = np.array([strategies.floats(min_value=0, max_value=box_length * 0.5,
                                          allow_nan=False, allow_infinity=False) for box_length in
                        periodic_box_lengths])
    point1 = np.array([draw(coord) for coord in lengths])
    point2 = np.array([draw(coord) for coord in lengths])

    # calculate the actual difference and magnitudes.
    diff = point1 - point2
    dist_sqr = np.dot(diff, diff)

    # move points to new random positions around the periodic box.
    images = strategies.integers(min_value=-100, max_value=100)
    images1 = np.array([draw(images) for x in range(3)])
    images2 = np.array([draw(images) for x in range(3)])
    point1_periodic = point1 + images1 * periodic_box_lengths
    point2_periodic = point2 + images2 * periodic_box_lengths

    return point1_periodic, point2_periodic, periodic_box_lengths, diff, dist_sqr


@given(points_in_periodic_box())
@example((np.zeros((3,)), np.zeros((3,)), np.array([1, 1, 1], dtype=float), np.zeros((3,)), 0))
@example((np.zeros((3,)), np.array([1, 1, 1]), np.array([2, 2, 2], dtype=float), np.array([-1, -1, -1]), 3))
@example((np.zeros((3,)), np.array([3, 3, 3]), np.array([2, 2, 2], dtype=float), np.array([1, 1, 1]), 3))
@example((np.zeros((3,)), np.array([0, 0.5, 0]), np.array([1, 1, 1], dtype=float), np.array([0, 0.5, 0]), 0.25))
def test_distance(vec_pbc_diff):
    point1, point2, periodic_box_lengths, expected_diff, expected_dist_sqr = vec_pbc_diff
    diff, dist_sqr = _calculate_diff_and_sqr_distance(point1, point2, periodic_box_lengths)
    assert np.allclose(dist_sqr, expected_dist_sqr, atol=1.0e-10)
    assert np.allclose(np.abs(diff), np.abs(expected_diff), atol=1.0e-10)
