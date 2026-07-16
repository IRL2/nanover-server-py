import numpy as np
from hypothesis import strategies as st, given

from MDAnalysis.lib import transformations
from nanover.utilities.transforms import (
    Transform,
    find_transformation_between_points,
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


@given(transformation=transformation())
def test_cube_alignment_valid(transformation):
    transform = Transform.from_parent_to_local_matrix(transformation)

    a = CUBE_POINTS
    b = np.array([transform.point_parent_to_local(point) for point in a])
    guess = Transform.from_parent_to_local_matrix(
        find_transformation_between_points(a, b)
    )
    c = np.array([guess.point_parent_to_local(point) for point in a])

    assert np.allclose(b, c)
