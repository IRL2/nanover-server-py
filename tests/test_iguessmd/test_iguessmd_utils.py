"""
Tests for :mod:`nanover.iguessmd.utils`.

Things to test for utility functions:
- get_every_nth returns arrays of the expected shapes with the expected values [√]
- calculate_unit_tangent correctly calculates the unit tangents for a given non-linear
  path
- calculate_com correctly calculates the COMs for a set of positions with known COMs
"""

import pytest
import numpy as np

from nanover.iguessmd.utils import get_every_nth, calculate_unit_tangents, calculate_com

from .iguessmd_test_utilities import define_circular_path

# Test systems for COM calculations, formatted as (positions, masses, expected COM)
TEST_COM_TWO_ATOMS = (
    np.array([[-1.0, 0.0, 0.0], [1.0, 0.0, 0.0]]),
    np.array([1.0, 1.0]),
    np.array([0.0, 0.0, 0.0]),
)
TEST_COM_METHANE = (
    np.array(
        [
            [1.0, 1.0, 0.0],
            [0.0, 1.0, -1.0 / np.sqrt(2)],
            [2.0, 1.0, -1.0 / np.sqrt(2)],
            [1.0, 0.0, 1.0 / np.sqrt(2)],
            [1.0, 2.0, 1.0 / np.sqrt(2)],
        ]
    ),
    np.array([12.01, 1.00, 1.00, 1.00, 1.00]),
    np.array([1.0, 1.0, 0.0]),
)
TEST_COM_CIRCLE = (
    np.array(
        [[np.cos(i * (np.pi / 4)), 0.0, np.sin(i * (np.pi / 4))] for i in range(8)]
    ),
    np.array([2.0, 1.05, 2.0, 1.05, 2.0, 1.05, 2.0, 1.05]),
    np.array([0.0, 0.0, 0.0]),
)
TEST_COM_CUBE = (
    np.array(
        [
            [0.0, 1.0, 2.0],
            [0.0, 1.0, 4.0],
            [0.0, 3.0, 2.0],
            [0.0, 3.0, 4.0],
            [2.0, 1.0, 2.0],
            [2.0, 1.0, 4.0],
            [2.0, 3.0, 2.0],
            [2.0, 3.0, 4.0],
        ]
    ),
    np.array([6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0]),
    np.array([1.0, 2.0, 3.0]),
)
TEST_iGUESSMD_CIRCULAR_PATH, TEST_iGUESSMD_CIRCULAR_PATH_TANGENTS = define_circular_path(10)

@pytest.mark.parametrize(
    "array,axis",
    [
        (np.arange(0, 99, 1), 0),
        (np.arange(0, 100, 1), 0),
        (np.arange(0, 101, 1), 0),
        (np.random.rand(3, 150, 4), 0),
        (np.random.rand(3, 150, 4), 1),
        (np.random.rand(3, 150, 4), 2),
    ],
)
@pytest.mark.parametrize(
    "every_nth, should_raise",
    [
        (-1, True),
        (0, True),
        (1, False),
        (2, False),
        (5, False),
        (10, False),
        (23, False),
    ],
)
def test_get_every_nth(array, axis, every_nth, should_raise):
    """
    Test that the function get_every_nth behaves as expected.
    """
    # Check whether invalid integer values raise ValueError
    if should_raise:
        with pytest.raises(ValueError):
            get_every_nth(array, axis=axis, every_nth=every_nth, include_end=False)
            assert every_nth <= 0
        return

    # Calculate reduced arrays
    reduced_array_without_end = get_every_nth(
        array, axis=axis, every_nth=every_nth, include_end=False
    )
    reduced_array_with_end = get_every_nth(
        array, axis=axis, every_nth=every_nth, include_end=True
    )

    # Calculate expected shapes
    expected_shape_without_end = np.round(
        np.floor((array.shape[axis] - 1) / every_nth) + 1
    )
    if (array.shape[axis] - 1) % every_nth != 0:
        expected_shape_with_end = expected_shape_without_end + 1
    else:
        expected_shape_with_end = expected_shape_without_end

    # Check dimensions match expected
    assert reduced_array_without_end.shape[axis] == expected_shape_without_end
    assert reduced_array_with_end.shape[axis] == expected_shape_with_end

    # Check final entries with and without end point included
    assert (
        np.take(reduced_array_with_end, -1, axis=axis) == np.take(array, -1, axis=axis)
    ).all()
    if (array.shape[axis] - 1) % every_nth != 0:
        assert (
            np.take(reduced_array_without_end, -1, axis=axis)
            != np.take(array, -1, axis=axis)
        ).all()


def test_calculate_unit_tangents_circular_path():
    """
    Test that checks whether the normalised tangent vectors of a circular
    3D path are calculated correctly, using the forward difference
    approximation and with the final tangent vector equal to the
    penultimate tangent vector.
    """

    calculated_tangent_vectors = calculate_unit_tangents(TEST_iGUESSMD_CIRCULAR_PATH)

    # Check that calculated tangents are as expected
    assert calculated_tangent_vectors.shape == TEST_iGUESSMD_CIRCULAR_PATH_TANGENTS.shape
    assert np.allclose(
        calculated_tangent_vectors, TEST_iGUESSMD_CIRCULAR_PATH_TANGENTS, rtol=1e-8
    )

    # Check calculated tangent vectors obey symmetry of circle (excluding final tangent)
    n_points = TEST_iGUESSMD_CIRCULAR_PATH.shape[0]
    assert np.allclose(
        calculated_tangent_vectors[: int((n_points - 1) / 2)],
        -calculated_tangent_vectors[int((n_points - 1) / 2): -1],
        rtol=1e-8,
    )


@pytest.mark.parametrize(
    "positions, masses, com",
    [TEST_COM_TWO_ATOMS, TEST_COM_METHANE, TEST_COM_CIRCLE, TEST_COM_CUBE],
)
def test_calculate_com(positions, masses, com):
    """
    Check that the function calculate_com correctly calculates
    the centre of mass of a set of atoms, given their positions
    and masses.
    """
    calculated_com = calculate_com(positions, masses)
    expected_com = com
    assert np.allclose(calculated_com, expected_com, atol=1e-16)
