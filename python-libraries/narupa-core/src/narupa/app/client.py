"""
Module containing a basic interactive moleculary dynamics client that receives frames
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

    def __init__(self, address:Optional[str]=None,
                 trajectory_port:Optional[int]=None,
                 imd_port:Optional[int]=None,
                 max_frames=50,
                 run_imd=True):
        self._frame_client = FrameClient(address=address, port=trajectory_port)
        self.max_frames = max_frames
        if run_imd:
            self._imd_client = ImdClient(address=address, port=imd_port)
        else:
            self._imd_client = None

        self._frames = deque(maxlen=self.max_frames)
        self._join_trajectory()

    @property
    def running_imd(self) -> bool:
        return self._imd_client is not None

    @property
    def latest_frame(self) -> FrameData:
        return self.frames[-1]

    @property
    def frames(self) -> Sequence[FrameData]:
        return self._frames

    def start_interaction(self) -> int:
        if self._imd_client is None:
            raise ValueError("Client started without IMD, cannot start an interaction!")
        return self._imd_client.start_interaction()

    def update_interaction(self, interaction_id, interaction):
        if self._imd_client is None:
            raise ValueError("Client started without IMD, cannot update an interaction!")
        self._imd_client.update_interaction(interaction_id, interaction)

    def stop_interaction(self, interaction_id) -> InteractionEndReply:
        if self._imd_client is None:
            raise ValueError("Client started without IMD, cannot stop an interaction!")
        return self._imd_client.stop_interaction(interaction_id)

    def _join_trajectory(self):
        self._frame_client.subscribe_frames_async(self._on_frame_received)

    def _on_frame_received(self, frame_index:int, frame:GrpcFrameData):
        self._frames.append(FrameData(frame))






