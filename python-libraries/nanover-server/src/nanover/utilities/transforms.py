import numpy as np
import numpy.typing as npt

from MDAnalysis.lib import transformations


class Transform:
    @classmethod
    def from_scene_pose(cls, pose: list[float]):
        tx, ty, tz, rx, ry, rz, rw, sx, sy, sz = pose

        translation = transformations.translation_matrix((tx, ty, tz))
        rotation = transformations.quaternion_matrix((rw, rx, ry, rz))
        scale = np.diagflat((-sx, sy, sz, 1.0))

        # compose in TRS order
        return cls.from_local_to_parent_matrix(translation @ rotation @ scale)

    @classmethod
    def from_local_to_parent_matrix(cls, local_to_parent: npt.NDArray):
        return cls(
            local_to_parent=local_to_parent,
            parent_to_local=np.linalg.inv(local_to_parent),
        )

    @classmethod
    def from_parent_to_local_matrix(cls, parent_to_local: npt.NDArray):
        return cls(
            parent_to_local=parent_to_local,
            local_to_parent=np.linalg.inv(parent_to_local),
        )

    def __init__(
        self,
        *,
        local_to_parent: npt.NDArray,
        parent_to_local: npt.NDArray,
    ):
        self._local_to_parent = local_to_parent
        self._parent_to_local = parent_to_local

    def point_local_to_parent(self, point):
        return _transform_vec3(self._local_to_parent, point)

    def point_parent_to_local(self, point):
        return _transform_vec3(self._parent_to_local, point)


def _transform_vec3(matrix, vector):
    return (matrix @ np.array([*vector[:3], 1]).reshape(4, 1)).reshape(4)[:3]


def find_transformation_between_point_patterns(
    start_points: npt.NDArray,
    final_points: npt.NDArray,
):
    """
    Return the 4x4 transformation matrix that best maps the start points to the final points.
    """
    # https://web.stanford.edu/class/cs273/refs/umeyama.pdf
    # https://doi.org/10.1109/34.88573
    p = start_points.T
    q = final_points.T

    # centroids
    cen_P = np.mean(p, axis=1).reshape(-1, 1)
    cen_Q = np.mean(q, axis=1).reshape(-1, 1)

    # centered vectors
    X = p - cen_P
    Y = q - cen_Q

    # svd
    U, sigma, Vt = np.linalg.svd(X @ Y.T)

    # reflection correction
    d = np.identity(Vt.T.shape[1])
    d[-1, -1] = np.linalg.det(Vt.T @ U.T)

    # rotation and translation
    R = Vt.T @ d @ U.T
    t = cen_Q - R @ cen_P

    # convert to TRS matrix4x4
    transform = np.identity(4)
    transform[:3, 3] = t.reshape(-1)
    transform[:3, :3] = R

    return transform
