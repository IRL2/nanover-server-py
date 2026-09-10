from os import PathLike

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
    index_entries: list[RecordingIndexEntry] = []
    readers: list[NanoverTrajectory] = []
    first_particle_frame = FrameData()
    first_frame = last_frame = None

    def frame_begins_next_universe(frame: FrameData):
        return frame.frame_dict.get(keys.FRAME_INDEX, None) == 0

    def finalise_prev_universe():
        nonlocal first_particle_frame, first_frame, last_frame

        reader = NanoverTrajectory.from_components(
            path=path,
            first_frame=first_particle_frame,
            index=list(index_entries),
            name=f"{path}[{first_frame}:{last_frame}]",
        )
        readers.append(reader)

        index_entries.clear()
        first_particle_frame = FrameData()
        first_frame = last_frame = None

    with NanoverRecordingReader.from_path(path) as reader:
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

    return readers


class NanoverTrajectory:
    @classmethod
    def from_components(
        cls,
        *,
        path: str | PathLike[str],
        first_frame: FrameData,
        index: list[RecordingIndexEntry],
        name: str = "Unnamed",
    ):
        reader = NanoverRecordingReader.from_path(path)
        reader.index = index
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
            self.reader, filename=self.name, convert_units=convert_units
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

    def __getitem__(self, key: slice | int):
        if isinstance(key, int):
            entry = self.reader.index[key]
            frame = self.first_frame.copy()
            frame.update(self.reader.get_frame_from_entry(entry))
            return frame
        elif isinstance(key, slice):
            entries = self.reader.index[key]

            def iterate():
                for entry in entries:
                    frame = self.first_frame.copy()
                    frame.update(self.reader.get_frame_from_entry(entry))
                    yield frame

            return iterate()

    def __len__(self):
        return len(self.reader)

    def __iter__(self):
        current = FrameData()
        current.update(self.first_frame)

        for entry, frame in self.reader.iter_frame_updates():
            current.update(frame)
            yield current.copy()

    def __repr__(self):
        return f"<NanoverTrajectory {self.name} with {len(self.reader)} frames of {self.first_frame.particle_count} atoms>"


def is_valid_first_frame(frame: FrameData):
    return all(key in frame for key in FIRST_FRAME_REQUIRED)
