import numpy as np
from nanover.utilities.transforms import Transform

from . import Mode, NanoverJupyterUtilities


def use_transform_handles(utilities: NanoverJupyterUtilities):
    cursor_grabbed_handle: dict[str, dict] = {}
    cursor_grabbed_matrix: dict = {}
    cursor_grabbed_original: dict = {}

    def intersect_all_handles(point):
        for key, handle in utilities.handles.all_prefixed_items():
            object_to_root = utilities.transforms.fetch_transform_root(handle["parent"])
            local_point = object_to_root.points_parent_to_local(point)
            center, radius = handle["sphere"]
            if np.linalg.norm(np.subtract(local_point, center)) < radius:
                return handle
        return None

    def cursor_in_object_parent_matrix(cursor: dict, object: str):
        object_parent = utilities.transforms.get_parent(object, default="root")
        parent_to_root = utilities.transforms.fetch_transform_root(
            object_parent
        ).local_to_parent_matrix
        cursor_to_root = Transform.from_state_cursor(cursor).local_to_parent_matrix
        return np.linalg.inv(parent_to_root) @ cursor_to_root

    def filter_matrix(cursor, target, handle: dict):
        cursor = cursor.copy()
        if not handle.get("translate", False):
            cursor[:3, 3] = target[:3, 3]
        return cursor

    class MoveObjectMode(Mode):
        def on_button_pressed(self, *, key: str, cursor: dict, button: str):
            hovered = intersect_all_handles(cursor["position"])
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

                object_in_parent = filter_matrix(
                    object_in_parent, cursor_grabbed_original.get(key), grabbed
                )

                utilities.transforms.update_transform(
                    transform,
                    transform=Transform.from_local_to_parent_matrix(object_in_parent),
                    parent=object_parent,
                )

            # show/remove hover graphic is this cursor hovers something
            hovered = grabbed or intersect_all_handles(cursor["position"])
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
