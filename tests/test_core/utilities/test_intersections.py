import math
from dataclasses import dataclass

import numpy as np
import numpy.typing as npt
from hypothesis import given, strategies as st
from nanover.utilities.intersection import ClosestPointResult, closest_point_on_polyline
from nanover.utilities.transforms import Transform

from .test_transforms import transformation

POINT_COUNT = 32
ARC_POINTS = np.array(
    [(math.cos(i / 32 * math.pi), math.sin(i / 32 * math.pi), 0) for i in range(32)]
)


@st.composite
def arcs(draw):
    transform = Transform.from_parent_to_local_matrix(draw(transformation()))
    return [transform.point_parent_to_local(point) for point in ARC_POINTS]


@dataclass(kw_only=True)
class Case:
    arc: list[npt.NDArray]
    query: npt.NDArray
    answer: ClosestPointResult


@st.composite
def cases(draw):
    # randomly transformed arc
    arc = draw(arcs())

    # pick segment (index), percentage along (t), and distance
    index = draw(st.integers(min_value=0, max_value=len(arc) - 2))
    t = draw(st.floats(min_value=0, max_value=1, allow_nan=False, allow_infinity=False))
    distance = draw(
        st.floats(min_value=0, max_value=10, allow_nan=False, allow_infinity=False)
    )

    # desired closest point is factor t between a and b
    a, b = arc[index], arc[index + 1]
    point = np.add(a, np.subtract(b, a) * t)

    # find point distance from desired closest point
    perp = np.cross(np.subtract(arc[1], arc[0]), np.subtract(arc[1], arc[2]))
    out = np.cross(perp, np.subtract(b, a))
    out = out / np.linalg.norm(out)
    query = point + out * distance

    return Case(
        arc=arc,
        query=query,
        answer=ClosestPointResult(
            point=point,
            index=index,
            distance=distance,
            t=t,
        ),
    )


@given(test_case=cases())
def test_closest_point_on_polyline(test_case: Case):
    answer = closest_point_on_polyline(test_case.arc, test_case.query)
    assert np.allclose(answer.point, test_case.answer.point, atol=0.001)
    assert np.allclose(answer.distance, test_case.answer.distance, atol=0.001)
