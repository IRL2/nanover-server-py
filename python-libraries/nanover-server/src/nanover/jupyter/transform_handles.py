import numpy as np
from nanover.utilities.transforms import (
    Transform,
    matrix_from_state_transform,
    state_transform_from_matrix,
)

from . import Mode, NanoverJupyterUtilities


def use_transform_handles(utilities: NanoverJupyterUtilities):
    """
    Add a mode that highlights hovered handles and allows grabbing and manipulation of them to manipulate the
    corresponding transform by pressing the primary button.
    """
    cursor_grabbed_handle: dict[str, dict] = {}
    cursor_grabbed_matrix: dict = {}
    cursor_grabbed_original: dict = {}

    def cursor_in_object_parent_matrix(cursor: dict, object: str):
        """
        Return the matrix of the cursor relative to the parent of given transform.
        """

        object_parent = utilities.transforms.get_parent(object, default="root")
        parent_to_root = utilities.transforms.fetch_transform_root(
            object_parent
        ).local_to_parent_matrix
        cursor_to_root = Transform.from_state_cursor(cursor).local_to_parent_matrix
        return np.linalg.inv(parent_to_root) @ cursor_to_root

    def filter_matrix_update(next_matrix, prev_matrix, handle: dict):
        """
        Return the next matrix with any translation/rotation/scale disallowed by the handle reverted to those of the
        previous matrix.
        """

        prev_pose = state_transform_from_matrix(prev_matrix)
        next_pose = list(state_transform_from_matrix(next_matrix))

        if not handle.get("translate", False):
            next_pose[0:3] = prev_pose[0:3]
        if not handle.get("rotate", False):
            next_pose[3:7] = prev_pose[3:7]
        if not handle.get("scale", False):
            next_pose[7:10] = prev_pose[7:10]

        return matrix_from_state_transform(next_pose)

    class MoveObjectMode(Mode):
        def on_button_pressed(self, *, key: str, cursor: dict, button: str):
            hovered = utilities.intersect_transform_handles(cursor["position"])
            available = (
                hovered not in cursor_grabbed_handle.values() and hovered is not None
            )

            # grab hovered object if not already grabbed
            if button == "primary" and available:
                transform = hovered["parent"]

                # cursor matrix relative to object parent
                cursor_in_parent = cursor_in_object_parent_matrix(cursor, transform)
                # object matrix relative to object parent
                object_in_parent = utilities.transforms.fetch_transform(
                    transform
                ).local_to_parent_matrix
                # matrix transforming cursor to object
                offset_matrix = np.linalg.inv(cursor_in_parent) @ object_in_parent

                cursor_grabbed_handle[key] = hovered
                cursor_grabbed_matrix[key] = offset_matrix
                cursor_grabbed_original[key] = object_in_parent

        def on_button_released(self, *, key: str, cursor: dict, button: str):
            # release grabbed
            if button == "primary":
                cursor_grabbed_handle.pop(key, None)
                cursor_grabbed_matrix.pop(key, None)
                cursor_grabbed_original.pop(key, None)

        def on_cursor_updated(self, *, key: str, cursor: dict):
            # if this cursor has grabbed an object, update the object pose from cursor pose
            grabbed = cursor_grabbed_handle.get(key, None)
            if grabbed is not None:
                transform = grabbed["parent"]

                object_parent = utilities.transforms.get_parent(
                    transform, default="root"
                )
                # cursor matrix relative to object parent
                cursor_in_parent = cursor_in_object_parent_matrix(cursor, transform)
                # matrix transforming cursor to object
                offset_matrix = cursor_grabbed_matrix.get(key, np.identity(4))
                # object matrix relative to object parent
                object_in_parent = cursor_in_parent @ offset_matrix

                object_in_parent = filter_matrix_update(
                    object_in_parent, cursor_grabbed_original.get(key), grabbed
                )

                utilities.transforms.update_transform(
                    transform,
                    transform=Transform.from_local_to_parent_matrix(object_in_parent),
                    parent=object_parent,
                )

            # show/remove hover graphic is this cursor hovers something
            hovered = grabbed or utilities.intersect_transform_handles(
                cursor["position"]
            )
            if hovered is None:
                utilities.objects.update_shape(f"hovered.{key}")
            else:
                center, radius = hovered["sphere"]
                utilities.objects.update_shape(
                    f"hovered.{key}",
                    position=center,
                    color=[1.0, 1.0, 0.0, 0.5],
                    size=radius * 2,
                    parent=hovered["parent"],
                )

    utilities.add_interaction_mode(MoveObjectMode, "move object", icon="✊")
