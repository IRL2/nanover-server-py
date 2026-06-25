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
        return cls(local_to_parent=translation @ rotation @ scale)

    def __init__(
        self,
        *,
        local_to_parent: npt.NDArray[float] | None = None,
        parent_to_local: npt.NDArray[float] | None = None,
    ):
        assert local_to_parent is None or parent_to_local is None
        self._local_to_parent = (
            local_to_parent
            if local_to_parent is not None
            else np.linalg.inv(parent_to_local)
        )
        self._parent_to_local = (
            parent_to_local
            if parent_to_local is not None
            else np.linalg.inv(local_to_parent)
        )

    def point_local_to_parent(self, point):
        return _transform_vec3(self._local_to_parent, point)

    def point_parent_to_local(self, point):
        return _transform_vec3(self._parent_to_local, point)


def _transform_vec3(matrix, vector):
    return (matrix @ np.array([*vector[:3], 1]).reshape(4, 1)).reshape(4)[:3]
