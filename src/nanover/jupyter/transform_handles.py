import math
from dataclasses import dataclass

import numpy as np

from . import Mode, NanoverJupyterUtilities
from .transform_grabbing import TransformGrabbingContext

_HOVER_CIRCLE_POINTS = 16
_HOVER_CIRCLE_ANGLES = [
    math.pi * 2 * i / (_HOVER_CIRCLE_POINTS - 1) for i in range(_HOVER_CIRCLE_POINTS)
]

_HOVER_CIRCLE = np.array(
    [[math.cos(angle) / 2, 0, math.sin(angle) / 2] for angle in _HOVER_CIRCLE_ANGLES]
)


@dataclass(kw_only=True)
class HandleGrabData:
    handle: dict


def use_transform_handles(utilities: NanoverJupyterUtilities):
    """
    Add a mode that highlights hovered handles and allows grabbing and manipulation of them to manipulate the
    corresponding transform by pressing the primary button.
    """

    context = TransformGrabbingContext[HandleGrabData].from_utilities(utilities)

    class MoveObjectMode(Mode):
        def on_button_pressed(self, *, key: str, cursor: dict, button: str):
            hovered = utilities.intersect_transform_handles(cursor["position"])
            # try grab hovered handle
            if button == "primary" and hovered is not None:
                grab = context.start_grab_from_cursor(
                    key, transform_id=hovered["parent"], cursor=cursor
                )
                if grab is not None:
                    grab.translate = hovered.get("translate", True)
                    grab.rotate = hovered.get("rotate", True)
                    grab.scale = hovered.get("scale", True)
                    grab.data = HandleGrabData(handle=hovered)

        def on_button_released(self, *, key: str, cursor: dict, button: str):
            # release cursor's grab
            if button == "primary":
                context.end_grab_from_cursor(key, cursor=cursor)

        def on_cursor_updated(self, *, key: str, cursor: dict):
            grab = context.update_grab_from_cursor(key, cursor=cursor)

            # show/remove hover graphic is this cursor hovers something
            hovered = (
                grab.data.handle
                if grab is not None
                else utilities.intersect_transform_handles(cursor["position"])
            )
            if hovered is None:
                utilities.objects.remove_line(f"hovered.{key}")
                utilities.objects.remove_shape(f"hovered.{key}")
            else:
                center, radius = hovered["sphere"]

                if hovered.get("minimal", False):
                    utilities.objects.update_line(
                        f"hovered.{key}",
                        positions=_HOVER_CIRCLE * radius + center,
                        color=[1.0, 1.0, 0.0, 0.5],
                        parent=hovered["parent"],
                    )
                else:
                    utilities.objects.update_shape(
                        f"hovered.{key}",
                        position=center,
                        color=[1.0, 1.0, 0.0, 0.5],
                        size=radius * 2,
                        parent=hovered["parent"],
                    )

    utilities.modes.add_mode(MoveObjectMode(), "move object", icon="✊")
