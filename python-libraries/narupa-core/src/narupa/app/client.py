# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.
"""
Module containing a basic interactive molecular dynamics client that receives frames
and can publish interactions.
"""
import time
from collections import deque
from typing import Optional, Sequence, Dict

from narupa.imd.particle_interaction import ParticleInteraction
from narupa.protocol.imd import InteractionEndReply
from narupa.trajectory import FrameClient, FrameData
from narupa.imd import ImdClient
from narupa.multiplayer import MultiplayerClient
from google.protobuf.struct_pb2 import Value

# Default to a low framerate to avoid build up in the frame stream
DEFAULT_SUBSCRIPTION_INTERVAL = 1 / 30


class NarupaClient:
    """
    Basic interactive molecular dynamics client that receives frames and can publish interactions.
    :param address Address of the Narupa frame server.
    :param trajectory_port: Port at which the Narupa frame server is running.
    :param imd_port: Port at which the Narupa IMD server is running.
    :param multiplayer_port: Port at which the Narupa multiplayer server is running.
    :param max_frames: Maximum number of frames to retain. If more frames are received, older frames will be removed.
    :param run_imd: Whether to connect to the IMD server.
    :param run_multiplayer: Whether to connect to the multiplayer server.
    :param all_frames: Whether to receive all frames produced by the trajectory server, or just subscribe to the latest
    frame.
    """
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
        """
        Closes the connection with the server.
        :param clear_frames: Whether to clear the frames received by the client, or keep them.
        """
        if clear_frames:
            self._first_frame = None
            self._frames.clear()
        self._frame_client.close()
        if self.running_imd:
            self._imd_client.close()

    def connect_trajectory(self, address: str, port: Optional[int] = None) -> None:
        """
        Connects the client to the given trajectory server, and begin receiving frames.
        :param address: The URL or IP address of the trajectory server.
        :param port: The port of the trajectory server.
        """
        self._frame_client = FrameClient(address=address, port=port)
        self._join_trajectory()

    def connect_imd(self, address: str, port: Optional[int] = None) -> None:
        """
        Connects the client to the given interactive molecular dynamics server,
        allowing it to start publishing interactions.
        :param address: The address of the IMD server.
        :param port: The port of the IMD server.
        """
        self._imd_client = ImdClient(address=address, port=port)

    def connect_multiplayer(self, address: str, port: Optional[int] = None) -> None:
        """
        Connects the client to the given multiplayer server.
        :param address: The address of the multiplayer server.
        :param port: The port of the multiplayer server.
        """
        self._multiplayer_client = MultiplayerClient(address=address, port=port)

    def connect(self, *, address=None, trajectory_port=None, imd_port=None,
                run_imd=True, multiplayer_port=None, run_multiplayer=True):
        """
        Connects the client to all the services, assuming they are running at the same URL on different ports.
        :param address: The address of the server(s).
        :param trajectory_port: The frame server port.
        :param imd_port: The IMD server port.
        :param run_imd: Whether to connect to the IMD server.
        :param multiplayer_port: The multiplayer server port.
        :param run_multiplayer: Whether to connect to the multiplayer server.
        """
        self.connect_trajectory(address, trajectory_port)
        if run_imd:
            self.connect_imd(address, imd_port)
        if run_multiplayer:
            self.connect_multiplayer(address, multiplayer_port)

    def wait_until_first_frame(self, check_interval=0.01, timeout=1):
        """
        Wait until the first frame is received from the server.
        :param check_interval: Interval at which to check if a frame has been received.
        :param timeout: Timeout after which to stop waiting for a frame.
        :return: The first :class:`FrameData` received.
        :raises: :class: Exception if no frame is received.
        """
        endtime = 0 if timeout is None else time.monotonic() + timeout

        while self.first_frame is None:
            if 0 < endtime < time.monotonic():
                raise Exception("Timed out waiting for first frame.")
            time.sleep(check_interval)

        return self.first_frame

    @property
    def running_imd(self) -> bool:
        """
        Indicates whether this client is running an IMD client or not.
        :return: `True` if this client is running an IMD client, `False` otherwise.
        """
        return self._imd_client is not None

    @property
    def running_multiplayer(self) -> bool:
        """
        Indicates whether this client is running a multiplayer client or not.
        :return: `True` if this client is running an multiplayer client, `False` otherwise.
        """
        return self._multiplayer_client is not None

    @property
    def latest_frame(self) -> Optional[FrameData]:
        """
        The trajectory frame most recently received, if any.
        :return: :class:`FrameData`, or `None` if none has been received.
        """
        if len(self.frames) is 0:
            return None
        else:
            return self.frames[-1]

    @property
    def latest_multiplayer_values(self) -> Dict[str, Value]:
        """
        The latest state of the multiplayer shared key/value store.
        :return: Dictionary of the current state of multiplayer shared key/value store.
        """
        return dict(self._multiplayer_client.resources)

    @property
    def frames(self) -> Sequence[FrameData]:
        """
        The most recently received frames up to the storage limit specified by `max_frames`.
        :return: Sequence of frames.
        """
        return list(self._frames)

    @property
    def first_frame(self) -> Optional[FrameData]:
        """
        The first received trajectory frame, if any.
        :return: The first frame received by this trajectory, or `None`.
        """
        return self._first_frame

    def start_interaction(self, interaction: Optional[ParticleInteraction] = None) -> int:
        """
        Start an interaction with the IMD server.
        :param interaction: An optional :class: ParticleInteraction with which to begin.
        :return: The unique interaction ID of this interaction, which can be used to update
        the interaction with :method: NarupaClient.update_interaction .

        :raises: ValueError, if the there is no IMD connection available.
        """
        if self._imd_client is None:
            raise ValueError("Client started without IMD, cannot start an interaction!")
        interaction_id = self._imd_client.start_interaction()
        if interaction is not None:
            self.update_interaction(interaction_id, interaction)
        return interaction_id

    def update_interaction(self, interaction_id, interaction:ParticleInteraction):
        """
        Updates the interaction identified with the given interaction_id on the server with
        parameters from the given interaction.
        :param interaction_id: The unique interaction ID, created with :method: NarupaClient.start_interaction,
        that identifies the interaction to update.
        :param interaction: The :class: ParticleInteraction providing new parameters for the interaction.

        :raises: ValueError, if the there is no IMD connection available, or if invalid parameters
        are passed to the server.
        :raises: KeyError, if the given interaction ID does not exist.
        """
        if self._imd_client is None:
            raise ValueError("Client started without IMD, cannot update an interaction!")
        self._imd_client.update_interaction(interaction_id, interaction)

    def stop_interaction(self, interaction_id) -> InteractionEndReply:
        """
        Stops the interaction identified with the given interaction_id on the server.
        :param interaction_id: The unique interaction ID, created with :method: NarupaClient.start_interaction,
        that identifies the interaction to stop.

        :raises: ValueError, if the there is no IMD connection available, or if invalid parameters
        are passed to the server.
        :raises: KeyError, if the given interaction ID does not exist.
        """
        if self._imd_client is None:
            raise ValueError("Client started without IMD, cannot stop an interaction!")
        return self._imd_client.stop_interaction(interaction_id)

    def join_multiplayer(self, player_name):
        """
        Joins multiplayer with the given player name.
        :param player_name: The player name with which to be identified.
        """
        self._multiplayer_client.join_multiplayer(player_name)

    def set_shared_value(self, key, value) -> bool:
        """
        Attempts to set the given key/value pair on the multiplayer shared value store.
        :param key: The key that identifies the value to be stored.
        :param value: The new value to store.
        :return: `True` if successful, `False` otherwise.
        """
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
