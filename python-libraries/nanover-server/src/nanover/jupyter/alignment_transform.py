from .utilities import TransformsUtility
from .frame_listener import FrameListener
from nanover.trajectory import FrameData
from nanover.utilities.transforms import StructureAlignment


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
        if self.key is not None and self.alignment is not None:
            transform = self.alignment.calculate_transform_to_framedata(full_frame)
            self.transforms.update_transform(
                self.key, transform=transform, parent="simulation"
            )
