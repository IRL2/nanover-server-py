from collections.abc import Sequence
from os import PathLike
from typing import Self, overload

from nanover.recording import NanoverRecordingReader
from nanover.recording.reading import RecordingIndexEntry
from nanover.trajectory import FrameData, keys

FIRST_FRAME_REQUIRED = {
    keys.PARTICLE_POSITIONS,
    keys.PARTICLE_COUNT,
}


def trajectories_from_recording(path: str | PathLike[str]):
    """
    Decompose a NanoVer trajectory into an mdanalysis Universe for each session of simulation (determined
    by frame_index resets).
    """
    return trajectories_from_reader(NanoverRecordingReader.from_path(path))


def trajectories_from_reader(reader: NanoverRecordingReader):
    index_entries: list[RecordingIndexEntry] = []
    trajectories: list[NanoverTrajectory] = []
    first_particle_frame = FrameData()
    first_frame = last_frame = None

    def frame_begins_next_universe(frame: FrameData):
        return frame.frame_dict.get(keys.FRAME_INDEX, None) == 0

    def finalise_prev_universe():
        nonlocal first_particle_frame, first_frame, last_frame

        trajectory = NanoverTrajectory.from_components(
            reader=reader.with_index(list(index_entries)),
            first_frame=first_particle_frame,
            name=f"{reader.name}[{first_frame}:{last_frame}]",
        )
        trajectories.append(trajectory)

        index_entries.clear()
        first_particle_frame = FrameData()
        first_frame = last_frame = None

    for i, entry in enumerate(reader):
        frame = reader.get_frame_from_entry(entry)

        if frame is None:
            continue

        if first_frame is None:
            first_frame = i
        last_frame = i

        if frame_begins_next_universe(frame) and index_entries:
            finalise_prev_universe()
        # aggregate initial frames until there is position and topology information
        if not is_valid_first_frame(first_particle_frame):
            first_particle_frame.update(frame)
        index_entries.append(entry)

    if is_valid_first_frame(first_particle_frame) and index_entries:
        finalise_prev_universe()

    return trajectories


class NanoverTrajectory(Sequence[FrameData]):
    @classmethod
    def from_components(
        cls,
        *,
        reader: NanoverRecordingReader,
        first_frame: FrameData,
        name: str = "Unnamed",
    ):
        return cls(
            reader=reader,
            first_frame=first_frame,
            name=name,
        )

    def to_universe(self, *, convert_units=True):
        from nanover.mdanalysis import frame_data_to_mdanalysis
        from nanover.mdanalysis.universe import NanoverReaderBase

        universe = frame_data_to_mdanalysis(self.first_frame)
        universe.trajectory = NanoverReaderBase(
            self.reader,
            first_frame=self.first_frame,
            filename=self.name,
            convert_units=convert_units,
        )
        return universe

    def __init__(
        self,
        *,
        reader: NanoverRecordingReader,
        first_frame: FrameData,
        name: str,
    ):
        self.name = name
        self.reader = reader
        self.first_frame = first_frame

    @overload
    def __getitem__(self, key: int) -> FrameData: ...

    @overload
    def __getitem__(self, key: slice) -> Self: ...

    def __getitem__(self, key: slice | int):
        if isinstance(key, int):
            entry = self.reader.index[key]
            frame = self.first_frame.copy()
            frame.update(self.reader.get_frame_from_entry(entry))
            return frame
        elif isinstance(key, slice):
            return self.from_components(
                reader=self.reader.with_index_sliced(key),
                first_frame=self.first_frame,
                name=f"{self.name}[{key}]",
            )

        raise TypeError("trajectory indices must be integers or slices, not type")

    def __len__(self):
        return len(self.reader)

    def __repr__(self):
        return f"<{self.__class__.__name__} {self.name} with {len(self.reader)} frames of {self.first_frame.particle_count} atoms>"


def is_valid_first_frame(frame: FrameData):
    return all(key in frame for key in FIRST_FRAME_REQUIRED)
