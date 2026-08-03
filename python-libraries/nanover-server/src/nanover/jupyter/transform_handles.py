import numpy as np
from nanover.utilities.transforms import Transform

from . import Mode, NanoverJupyterUtilities


def use_transform_handles(utilities: NanoverJupyterUtilities):
    cursor_grabbed_object: dict[str, str] = {}
    cursor_grabbed_matrix: dict = {}

    def intersect_all_transforms(point):
        for key in utilities.transforms.all_prefixed():
            object_to_root = utilities.transforms.fetch_transform_root(key)
            local_point = object_to_root.points_parent_to_local(point)
            if np.linalg.norm(local_point) < 0.25:
                return key
        return None

    def cursor_in_object_parent_matrix(cursor: dict, object: str):
        object_parent = utilities.transforms.get_parent(object, default="root")
        cursor_to_root = Transform.from_state_cursor(cursor).local_to_parent_matrix
        parent_to_root = utilities.transforms.fetch_transform_root(
            object_parent
        ).local_to_parent_matrix
        return np.linalg.inv(parent_to_root) @ cursor_to_root

    class MoveObjectMode(Mode):
        def on_button_pressed(self, *, key: str, cursor: dict, button: str):
            hovered = intersect_all_transforms(cursor["position"])
            available = (
                hovered not in cursor_grabbed_object.values() and hovered is not None
            )

            # grab hovered object if not already grabbed
            if button == "primary" and available:
                # cursor matrix relative to object parent
                cursor_in_parent = cursor_in_object_parent_matrix(cursor, hovered)
                # object matrix relative to object parent
                object_in_parent = utilities.transforms.fetch_transform(
                    hovered
                ).local_to_parent_matrix
                # matrix transforming cursor to object
                offset_matrix = np.linalg.inv(cursor_in_parent) @ object_in_parent

                cursor_grabbed_object[key] = hovered
                cursor_grabbed_matrix[key] = offset_matrix

        def on_button_released(self, *, key: str, cursor: dict, button: str):
            # release grabbed
            if button == "primary":
                cursor_grabbed_object.pop(key, None)
                cursor_grabbed_matrix.pop(key, None)

        def on_cursor_updated(self, *, key: str, cursor: dict):
            # if this cursor has grabbed an object, update the object pose from cursor pose
            grabbed = cursor_grabbed_object.get(key, None)
            if grabbed is not None:
                object_parent = utilities.transforms.get_parent(grabbed, default="root")
                # cursor matrix relative to object parent
                cursor_in_parent = cursor_in_object_parent_matrix(cursor, grabbed)
                # matrix transforming cursor to object
                offset_matrix = cursor_grabbed_matrix.get(key, np.identity(4))
                # object matrix relative to object parent
                object_in_parent = cursor_in_parent @ offset_matrix

                utilities.transforms.update_transform(
                    grabbed,
                    transform=Transform.from_local_to_parent_matrix(object_in_parent),
                    parent=object_parent,
                )

            # show/remove hover graphic is this cursor hovers something
            hovered = intersect_all_transforms(cursor["position"])
            if hovered is None:
                utilities.objects.update_shape(f"hovered.{key}")
            else:
                utilities.objects.update_shape(
                    f"hovered.{key}",
                    color=[1.0, 1.0, 0.0, 0.5],
                    size=0.25,
                    parent=hovered,
                )

    utilities.add_interaction_mode(MoveObjectMode, "move object", icon="✊")
