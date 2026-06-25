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
