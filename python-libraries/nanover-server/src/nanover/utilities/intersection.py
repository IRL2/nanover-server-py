from dataclasses import dataclass
import numpy as np


@dataclass(kw_only=True)
class ClosestPointResult:
    point: np.ndarray
    distance: float
    index: int
    t: float


def closest_point_on_polyline(points: list, target):
    # each point to next delta
    p2p1 = np.diff(points, axis=0)
    # each point to next length squared
    p2p1dot = (p2p1 * p2p1).sum(1)
    # each point to target delta
    ptp1 = np.subtract(target, points)[:-1]
    # percent along each segment the target is
    t = (ptp1 * p2p1).sum(1) / p2p1dot
    # clamp between 0,1
    np.clip(t, 0, 1, out=t)

    # each segment closest point (p1 + t * p2p1)
    pc = np.add(points[:-1], t.reshape(-1, 1) * p2p1)
    # each closest point to target delta
    ptpc = np.subtract(target, pc)
    # each closest point to target length squared
    ptpcdot = (ptpc * ptpc).sum(1)

    # index of closest point to target minimum length squared
    index = np.argmin(ptpcdot)

    return ClosestPointResult(
        point=points[index],  # TODO: should be pc[index]?
        distance=float(np.sqrt(ptpcdot[index])),
        index=int(index),
        t=t,
    )
