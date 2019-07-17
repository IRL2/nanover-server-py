"""
Module containing a basic interactive molecular dynamics client that receives frames
and can publish interactions.
"""
from collections import deque
from typing import Optional, Sequence

from narupa.imd.imd_client import ImdClient
from narupa.protocol.imd import InteractionEndReply
from narupa.trajectory import FrameClient, FrameData
from narupa.protocol.trajectory import FrameData as GrpcFrameData


class NarupaClient:
    _imd_client: ImdClient
    _frame_client: FrameClient
    _frames: deque

    def __init__(self, address: Optional[str] = None,
                 trajectory_port: Optional[int] = None,
                 imd_port: Optional[int] = None,
                 max_frames=50,
                 run_imd=True,
                 all_frames=True):
        self.all_frames = all_frames
        self.connect(address, trajectory_port, imd_port, run_imd)

        self.max_frames = max_frames
        self._frames = deque(maxlen=self.max_frames)
        self._first_frame = None

    def close(self, clear_frames=True):
        if clear_frames:
            self._first_frame = None
            self._frames.clear()
        self._frame_client.close()
        if self.running_imd:
            self._imd_client.close()

    def connect(self, address=None, trajectory_port=None, imd_port=None, run_imd=True):
        self._frame_client = FrameClient(address=address, port=trajectory_port)
        if run_imd:
            self._imd_client = ImdClient(address=address, port=imd_port)
        else:
            self._imd_client = None
        self._join_trajectory()


    @property
    def running_imd(self) -> bool:
        return self._imd_client is not None

    @property
    def latest_frame(self) -> FrameData:
        if len(self.frames) is 0:
            return None
        else:
            return self.frames[-1]

    @property
    def frames(self) -> Sequence[FrameData]:
        return self._frames

    @property
    def first_frame(self) -> FrameData:
        return self._first_frame

    def start_interaction(self, interaction=None) -> int:
        if self._imd_client is None:
            raise ValueError("Client started without IMD, cannot start an interaction!")
        interaction_id = self._imd_client.start_interaction()
        if interaction is not None:
            self.update_interaction(interaction_id, interaction)
        return interaction_id

    def update_interaction(self, interaction_id, interaction):
        if self._imd_client is None:
            raise ValueError("Client started without IMD, cannot update an interaction!")
        self._imd_client.update_interaction(interaction_id, interaction)

    def stop_interaction(self, interaction_id) -> InteractionEndReply:
        if self._imd_client is None:
            raise ValueError("Client started without IMD, cannot stop an interaction!")
        return self._imd_client.stop_interaction(interaction_id)

    def _join_trajectory(self):
        if self.all_frames:
            self._frame_client.subscribe_frames_async(self._on_frame_received)
        else:
            self._frame_client.subscribe_last_frames_async(self._on_frame_received)

    def _on_frame_received(self, frame_index: int, frame: FrameData):
        if self._first_frame is None:
            self._first_frame = frame
        self._frames.append(frame)
