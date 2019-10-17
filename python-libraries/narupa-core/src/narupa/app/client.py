# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.
"""
Module containing a basic interactive molecular dynamics client that receives frames
and can publish interactions.
"""
import time
from collections import deque
from typing import Optional, Sequence, Dict

from narupa.protocol.imd import InteractionEndReply
from narupa.trajectory import FrameClient, FrameData
from narupa.imd import ImdClient
from narupa.multiplayer import MultiplayerClient
from google.protobuf.struct_pb2 import Value

# Default to a low framerate to avoid build up in the frame stream
DEFAULT_SUBSCRIPTION_INTERVAL = 1 / 30


class NarupaClient:
    _frame_client: FrameClient
    _imd_client: ImdClient
    _multiplayer_client: MultiplayerClient
    _frames: deque

    def __init__(self, *,
                 address: Optional[str] = None,
                 trajectory_port: Optional[int] = None,
                 imd_port: Optional[int] = None,
                 multiplayer_port: Optional[int] = None,
                 max_frames=50,
                 run_imd=True,
                 run_multiplayer=True,
                 all_frames=True):
        self.all_frames = all_frames
        self.max_frames = max_frames

        self._frame_client = None
        self._imd_client = None
        self._multiplayer_client = None

        self.connect(address=address,
                     trajectory_port=trajectory_port,
                     imd_port=imd_port,
                     multiplayer_port=multiplayer_port,
                     run_imd=run_imd,
                     run_multiplayer=run_multiplayer)

        self._frames = deque(maxlen=self.max_frames)
        self._first_frame = None

    def close(self, clear_frames=True):
        if clear_frames:
            self._first_frame = None
            self._frames.clear()
        self._frame_client.close()
        if self.running_imd:
            self._imd_client.close()

    def connect_trajectory(self, address: str, port: Optional[int] = None) -> None:
        self._frame_client = FrameClient(address=address, port=port)
        self._join_trajectory()

    def connect_imd(self, address: str, port: Optional[int] = None) -> None:
        self._imd_client = ImdClient(address=address, port=port)

    def connect_multiplayer(self, address: str, port: Optional[int] = None) -> None:
        self._multiplayer_client = MultiplayerClient(address=address, port=port)

    def connect(self, *, address=None, trajectory_port=None, imd_port=None,
                run_imd=True, multiplayer_port=None, run_multiplayer=True):
        self.connect_trajectory(address, trajectory_port)
        if run_imd:
            self.connect_imd(address, imd_port)
        if run_multiplayer:
            self.connect_multiplayer(address, multiplayer_port)

    def wait_until_first_frame(self, check_interval=0.01, timeout=1):
        endtime = 0 if timeout is None else time.monotonic() + timeout

        while self.first_frame is None:
            if 0 < endtime < time.monotonic():
                raise Exception("Timed out waiting for first frame.")
            time.sleep(check_interval)

        return self.first_frame

    @property
    def running_imd(self) -> bool:
        return self._imd_client is not None

    @property
    def running_multiplayer(self) -> bool:
        return self._multiplayer_client is not None

    @property
    def latest_frame(self) -> Optional[FrameData]:
        """
        The trajectory frame most recently received, if any.
        """
        if len(self.frames) is 0:
            return None
        else:
            return self.frames[-1]

    @property
    def latest_multiplayer_values(self) -> Dict[str, Value]:
        """
        The latest state of the multiplayer shared key/value store.
        """
        return dict(self._multiplayer_client.resources)

    @property
    def frames(self) -> Sequence[FrameData]:
        """
        The most recently received frames up to some storage limit.
        """
        return self._frames

    @property
    def first_frame(self) -> Optional[FrameData]:
        """
        The first received trajectory, if any.
        """
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

    def join_multiplayer(self, player_name):
        self._multiplayer_client.join_multiplayer(player_name)

    def set_shared_value(self, key, value) -> bool:
        return self._multiplayer_client.try_set_resource_value(key, value)

    def _join_trajectory(self):
        if self.all_frames:
            self._frame_client.subscribe_frames_async(self._on_frame_received)
        else:
            self._frame_client.subscribe_last_frames_async(
                self._on_frame_received,
                DEFAULT_SUBSCRIPTION_INTERVAL,
            )

    def _on_frame_received(self, frame_index: int, frame: FrameData):
        if self._first_frame is None:
            self._first_frame = frame
        self._frames.append(frame)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
