import numpy as np
from nanover.trajectory import FrameData, MissingDataError
from nanover.utilities.transforms import StructureAlignment, Transform

from .frame_listener import FrameListener
from .utilities import TransformsUtility


class AlignmentTransform(FrameListener):
    key: str | None = None
    alignment: StructureAlignment | None = None

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.transforms = TransformsUtility(self._app_server)

    def config(self, *, key: str, alignment: StructureAlignment):
        self.key = key
        self.alignment = alignment

    def on_frame_update(self, full_frame: FrameData, frame_update: FrameData):
        if self.key is None or self.alignment is None:
            return

        try:
            transform = self.alignment.calculate_transform_to_framedata(full_frame)
        except MissingDataError:
            transform = Transform.from_parent_to_local_matrix(np.identity(4))

        self.transforms.update_transform(
            self.key, transform=transform, parent="simulation"
        )
