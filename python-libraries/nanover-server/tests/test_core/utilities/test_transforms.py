import numpy as np
from hypothesis import given, strategies as st
from MDAnalysis.lib import transformations
from nanover.utilities.transforms import (
    Transform,
    find_transformation_between_point_patterns,
)


def coords():
    return st.floats(
        allow_nan=False,
        allow_infinity=False,
        min_value=-1.0e8,
        max_value=1.0e8,
    )


def vec3s():
    return st.lists(coords(), min_size=3, max_size=3)


CUBE_POINTS = np.array(
    [(x, y, z) for x in range(2) for y in range(2) for z in range(2)]
)


@st.composite
def transformation(draw):
    translation = draw(vec3s())
    axis = draw(vec3s().filter(lambda vec: np.sum(vec) > 0.01))
    angle = draw(coords())

    translation = transformations.translation_matrix(translation)
    rotation = transformations.rotation_matrix(angle, axis)

    return translation @ rotation


@given(transformation=transformation(), points=st.lists(vec3s(), min_size=1))
def test_transform_points_equals_transform_point(transformation, points):
    transform = Transform.from_parent_to_local_matrix(transformation)

    a = transform.points_parent_to_local(points)
    b = [transform.point_parent_to_local(point) for point in points]

    assert np.allclose(a, b)

    a = transform.points_local_to_parent(points)
    b = [transform.point_local_to_parent(point) for point in points]

    assert np.allclose(a, b)


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
