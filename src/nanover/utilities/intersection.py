from dataclasses import dataclass

import numpy as np


@dataclass(kw_only=True)
class PointAlongResult:
    point: np.ndarray
    index: int
    t: float


def point_at_distance_along_polyline(points: list, distance):
    # each point to next delta
    p1p2 = np.diff(points, axis=0)
    # each point to next length
    p1p2l = np.linalg.norm(p1p2, axis=1)
    # each point distance along polyline
    p1p2d = np.cumsum(p1p2l)

    # None if beyond total line length
    if p1p2d[-1] < distance:
        return None

    # last point under distance
    index = np.searchsorted(p1p2d, distance)
    # distance after last point
    extra = distance - (p1p2d[index - 1] if index > 0 else 0)
    # percent along segment after last point
    t = extra / p1p2l[index]
    # point along segment after last point (p1 + t * p1p2)
    p = np.add(points[index], t * p1p2[index])
    return PointAlongResult(
        point=p,
        index=int(index),
        t=float(t),
    )


@dataclass(kw_only=True)
class ClosestPointResult:
    point: np.ndarray
    distance: float
    index: int
    t: float


def closest_point_on_polyline(points: list, target):
    # each point to next delta
    p1p2 = np.diff(points, axis=0)
    # each point to next length squared
    p1p2dot = (p1p2 * p1p2).sum(1)
    # each point to target delta
    ptp1 = np.subtract(target, points)[:-1]
    # percent along each segment the target is
    t = (ptp1 * p1p2).sum(1) / p1p2dot
    # clamp between 0,1
    np.clip(t, 0, 1, out=t)

    # each segment closest point (p1 + t * p1p2)
    pc = np.add(points[:-1], t.reshape(-1, 1) * p1p2)
    # each closest point to target delta
    ptpc = np.subtract(target, pc)
    # each closest point to target length squared
    ptpcdot = (ptpc * ptpc).sum(1)

    # index of closest point to target minimum length squared
    index = np.argmin(ptpcdot)

    return ClosestPointResult(
        point=pc[index],
        distance=float(np.sqrt(ptpcdot[index])),
        index=int(index),
        t=t[index],
    )
