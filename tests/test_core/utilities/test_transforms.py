import numpy as np
from hypothesis import given, strategies as st
from nanover.testing.strategies import transformation, vec3s
from nanover.utilities.transforms import (
    STATE_TRANSFORM_IDENTITY,
    Transform,
    find_transformation_between_point_patterns,
    matrix_from_state_transform,
    state_transform_from_matrix,
    unpack_partial_state_transform,
)

CUBE_POINTS = np.array(
    [(x, y, z) for x in range(2) for y in range(2) for z in range(2)]
)


@given(transformation=transformation(), points=st.lists(vec3s(), min_size=1))
def test_transform_points_equals_transform_point(transformation, points):
    transform = Transform.from_parent_to_local_matrix(transformation)

    a = transform.points_parent_to_local(points)
    b = [transform.point_parent_to_local(point) for point in points]

    assert np.allclose(a, b, atol=1e-7)

    a = transform.points_local_to_parent(points)
    b = [transform.point_local_to_parent(point) for point in points]

    assert np.allclose(a, b, atol=1e-7)


@given(transformation=transformation())
def test_cube_alignment_valid(transformation):
    transform = Transform.from_parent_to_local_matrix(transformation)

    a = CUBE_POINTS
    b = transform.points_parent_to_local(a)
    guess = Transform.from_parent_to_local_matrix(
        find_transformation_between_point_patterns(a, b)
    )
    c = guess.points_parent_to_local(a)

    assert np.allclose(b, c)


@given(transformation=transformation())
def test_unpack_partial(transformation):
    state_transform = state_transform_from_matrix(transformation)

    for i in range(10):
        unpacked = list(unpack_partial_state_transform(state_transform[:i]))
        expected = [*state_transform[:i], *STATE_TRANSFORM_IDENTITY[i:]]
        assert unpacked == expected


@given(transformation=transformation())
def test_transform_matrix_reflexive(transformation):
    converted = matrix_from_state_transform(state_transform_from_matrix(transformation))
    assert np.allclose(converted, transformation)
