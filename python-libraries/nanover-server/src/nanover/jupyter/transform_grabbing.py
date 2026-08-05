from dataclasses import dataclass

import numpy as np
import numpy.typing as npt
from nanover.utilities.transforms import (
    Transform,
    matrix_from_state_transform,
    state_transform_from_matrix,
)

from . import NanoverJupyterUtilities


@dataclass(kw_only=True)
class TransformGrabData[TData: object]:
    transform_id: str
    transform_initial_matrix: npt.NDArray
    offset_matrix: npt.NDArray

    translate: bool = True
    rotate: bool = True
    scale: bool = True

    data: TData | None


class TransformGrabbingContext[TData]:
    @classmethod
    def from_utilities(cls, utilities: NanoverJupyterUtilities):
        return cls(utilities)

    def __init__(self, utilities: NanoverJupyterUtilities):
        self.utilities = utilities
        self._grabs: dict[str, TransformGrabData[TData]] = {}

    def add_grab(self, key: str, grab: TransformGrabData[TData]):
        self._grabs[key] = grab

    def get_grab(self, key: str):
        return self._grabs.get(key, None)

    def end_grab(self, key: str):
        return self._grabs.pop(key, None)

    def start_grab_from_cursor(
        self,
        key: str,
        *,
        transform_id: str,
        cursor: dict,
        translate=True,
        rotate=True,
        scale=True,
    ):
        available = not any(
            transform_id == grab.transform_id for grab in self._grabs.values()
        )
        if not available:
            return None

        # cursor matrix relative to object parent
        cursor_in_parent = self.cursor_in_object_parent_matrix(cursor, transform_id)
        # object matrix relative to object parent
        object_in_parent = self.utilities.transforms.fetch_transform(
            transform_id
        ).local_to_parent_matrix
        # matrix transforming cursor to object
        offset_matrix = np.linalg.inv(cursor_in_parent) @ object_in_parent

        grab = TransformGrabData(
            transform_id=transform_id,
            transform_initial_matrix=object_in_parent,
            offset_matrix=offset_matrix,
            translate=translate,
            rotate=rotate,
            scale=scale,
            data=None,
        )

        self.add_grab(key, grab)

        return grab

    def end_grab_from_cursor(self, key: str, *, cursor: dict):
        return self.end_grab(key)

    def update_grab_from_cursor(self, key: str, *, cursor: dict):
        grab = self._grabs.get(key, None)
        if grab is None:
            return None

        object_parent = self.utilities.transforms.get_parent(
            grab.transform_id, default="root"
        )
        # cursor matrix relative to object parent
        cursor_in_parent = self.cursor_in_object_parent_matrix(
            cursor, grab.transform_id
        )
        # object matrix relative to object parent
        object_in_parent = cursor_in_parent @ grab.offset_matrix

        object_in_parent = constrain_matrix_update(
            object_in_parent,
            grab.transform_initial_matrix,
            translate=grab.translate,
            rotate=grab.rotate,
            scale=grab.scale,
        )

        self.utilities.transforms.update_transform(
            grab.transform_id,
            transform=Transform.from_local_to_parent_matrix(object_in_parent),
            parent=object_parent,
        )

        return grab

    def transform_is_grabbed(self, key):
        return any(key == grab.transform_id for grab in self._grabs.values())

    def cursor_in_object_parent_matrix(
        self,
        cursor: dict,
        object: str,
    ):
        """
        Return the matrix of the cursor relative to the parent of given transform.
        """

        object_parent = self.utilities.transforms.get_parent(object, default="root")
        root_to_parent = self.utilities.transforms.fetch_transform_root(
            object_parent
        ).parent_to_local_matrix
        cursor_to_root = Transform.from_state_cursor(cursor).local_to_parent_matrix
        return root_to_parent @ cursor_to_root


def constrain_matrix_update(
    next_matrix, prev_matrix, *, translate=True, rotate=True, scale=True
):
    """
    Return the next matrix with any translation/rotation/scale disallowed reverted to those of the previous matrix.
    """

    prev_pose = state_transform_from_matrix(prev_matrix)
    next_pose = list(state_transform_from_matrix(next_matrix))

    if not translate:
        next_pose[0:3] = prev_pose[0:3]
    if not rotate:
        next_pose[3:7] = prev_pose[3:7]
    if not scale:
        next_pose[7:10] = prev_pose[7:10]

    return matrix_from_state_transform(next_pose)
