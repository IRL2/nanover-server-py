# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.
"""
Module containing a basic interactive molecular dynamics client that receives frames
and can publish interactions.
"""
import time
from collections import deque, ChainMap
from functools import wraps
from typing import Iterable
from typing import Optional, Sequence, Dict, MutableMapping

from google.protobuf.struct_pb2 import Value
from grpc import RpcError, StatusCode
from narupa.app.selection import RenderingSelection
from narupa.core import CommandInfo, NarupaClient
from narupa.core.protobuf_utilities import struct_to_dict, dict_to_struct
from narupa.imd import ImdClient
from narupa.imd.particle_interaction import ParticleInteraction
from narupa.multiplayer import MultiplayerClient
from narupa.protocol.imd import InteractionEndReply
from narupa.trajectory import FrameClient, FrameData
from narupa.trajectory.frame_server import PLAY_COMMAND_KEY, STEP_COMMAND_KEY, PAUSE_COMMAND_KEY, RESET_COMMAND_KEY

# Default to a low framerate to avoid build up in the frame stream
DEFAULT_SUBSCRIPTION_INTERVAL = 1 / 30

# ID of the root selection
SELECTION_ROOT_ID = 'selection.root'
# Name of the root selection
SELECTION_ROOT_NAME = 'Root Selection'


def _update_commands(client: NarupaClient):
    try:
        return client.update_available_commands()
    except RpcError as e:
        if e._state.code == StatusCode.UNAVAILABLE:
            return {}
        else:
            raise e


class NarupaImdClient:
    """
    Interactive molecular dynamics client that receives frames, create selections,
    and join the multiplayer shared state.

    :param address Address of the Narupa frame server.
    :param trajectory_port: Port at which the Narupa frame server is running.
    :param imd_port: Port at which the Narupa IMD server is running.
    :param multiplayer_port: Port at which the Narupa multiplayer server is running.
    :param max_frames: Maximum number of frames to retain. If more frames are
        received, older frames will be removed.
    :param all_frames: Whether to receive all frames produced by the trajectory
        server, or just subscribe to the latest frame.


    Inspecting a Frame
    ==================

    The Narupa Imd client can be used to inspect frames received from a :class:`narupa.trajectory.FrameServer`,
    which can be useful for analysis.

    .. python
        # Assuming a server has been created with default address and ports.
        client = NarupaImdClient()
        # Fetch the first frame.
        first_frame = client.wait_until_first_frame(check_interval=0.5, timeout=10)
        # Print the number of particles in the frame
        print(first_frame.particle_count)

    Creating Selections for Rendering
    =================================

    One of the main uses of the client is to create selections that control how a group of particles
    will be visualised and interacted with in other clients (e.g., the VR client):

    .. code-block::  python

        # Connect to the multiplayer
        client.join_multiplayer("Selection Example")
        # Create a selection called 'Selection' which selects particles with indices 0-4
        selection = client.create_selection("Selection", [0, 1, 2, 3, 4])

    Selections are created and updated based on lists of particle indices. Tools such as
    `MDAnalysis <https://www.mdanalysis.org/>`_ or `MDTraj <http://mdtraj.org/>_` are very good at
    extracting indices of particles based on a human readable command.

    With a selection in hand, the way in which it is rendered and interacted with can be changed and
    transmitted to other clients (i.e. VR) using the `modify` context:

    .. code-block::  python

        # Change how the selection is rendered and interacted with.
        with selection.modify():
            selection.renderer = {
                    'color': 'IndianRed',
                    'scale': 0.1,
                    'render': 'liquorice'
                }
            selection.velocity_reset = True  # Reset the velocities after interacting.
            selection.interaction_method = 'group'  # Interact with the selection as a group.

    """
    _frame_client: FrameClient
    _imd_client: ImdClient
    _multiplayer_client: MultiplayerClient
    _frames: deque

    _next_selection_id: int = 0

    def __init__(self, *,
                 address: Optional[str] = None,
                 trajectory_port: Optional[int] = None,
                 imd_port: Optional[int] = None,
                 multiplayer_port: Optional[int] = None,
                 max_frames=50,
                 all_frames=True):
        self.all_frames = all_frames
        self.max_frames = max_frames

        self.connect(address=address,
                     trajectory_port=trajectory_port,
                     imd_port=imd_port,
                     multiplayer_port=multiplayer_port)

        self._frames = deque(maxlen=self.max_frames)
        self._first_frame = None

        self.update_available_commands()  # initialise the set of available commands.

    def close(self, clear_frames=True):
        """
        Closes the connection with the server.

        :param clear_frames: Whether to clear the frames received by the
            client, or keep them.
        """
        if clear_frames:
            self._first_frame = None
            self._frames.clear()
        self._frame_client.close()
        self._imd_client.close()
        self._multiplayer_client.close()

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
        self._join_interactions()

    def connect_multiplayer(self, address: str, port: Optional[int] = None) -> None:
        """
        Connects the client to the given multiplayer server.

        :param address: The address of the multiplayer server.
        :param port: The port of the multiplayer server.
        """
        self._multiplayer_client = MultiplayerClient(address=address, port=port)

    def connect(self, *, address=None, trajectory_port=None, imd_port=None,
                multiplayer_port=None):
        """
        Connects the client to all the services, assuming they are running at
        the same URL on different ports.

        :param address: The address of the server(s).
        :param trajectory_port: The frame server port.
        :param imd_port: The IMD server port.
        :param multiplayer_port: The multiplayer server port.
        """
        self.connect_trajectory(address, trajectory_port)
        self.connect_imd(address, imd_port)
        self.connect_multiplayer(address, multiplayer_port)

    def wait_until_first_frame(self, check_interval=0.01, timeout=1):
        """
        Wait until the first frame is received from the server.

        :param check_interval: Interval at which to check if a frame has been
            received.
        :param timeout: Timeout after which to stop waiting for a frame.
        :return: The first :class:`FrameData` received.
        :raises Exception: if no frame is received.
        """
        endtime = 0 if timeout is None else time.monotonic() + timeout

        while self.first_frame is None:
            if 0 < endtime < time.monotonic():
                raise Exception("Timed out waiting for first frame.")
            time.sleep(check_interval)

        return self.first_frame

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
        The most recently received frames up to the storage limit specified
        by `max_frames`.

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

    @property
    def interactions(self) -> Dict[str, ParticleInteraction]:
        """
        The dictionary of current interactions received by this client.
        :return: Dictionary of active interactions, keyed by interaction ID identifying who is performing the
        interactions.
        """
        return self._imd_client.interactions

    def start_interaction(self, interaction: Optional[ParticleInteraction] = None) -> str:
        """
        Start an interaction with the IMD server.

        :param interaction: An optional :class: ParticleInteraction with which
            to begin.
        :return: The unique interaction ID of this interaction, which can be
            used to update the interaction with
            :func:`~NarupaClient.update_interaction`.

        :raises: ValueError, if the there is no IMD connection available.
        """
        interaction_id = self._imd_client.start_interaction()
        if interaction is not None:
            self.update_interaction(interaction_id, interaction)
        return interaction_id

    def update_interaction(self, interaction_id, interaction: ParticleInteraction):
        """
        Updates the interaction identified with the given interaction_id on
        the server with parameters from the given interaction.

        :param interaction_id: The unique interaction ID, created with
            :func:`~NarupaClient.start_interaction`, that identifies the
            interaction to update.
        :param interaction: The :class: ParticleInteraction providing new
            parameters for the interaction.

        :raises: ValueError, if the there is no IMD connection available, or
            if invalid parameters are passed to the server.
        :raises KeyError: if the given interaction ID does not exist.
        """
        self._imd_client.update_interaction(interaction_id, interaction)

    def stop_interaction(self, interaction_id) -> InteractionEndReply:
        """
        Stops the interaction identified with the given interaction_id on the server.
        This method blocks until the server returns a reply indicating that the
        interaction has stopped.

        :param interaction_id: The unique interaction ID, created with
            :func:`~NarupaClient.start_interaction`, that identifies the
            interaction to stop.
        :return: An :class:`InteractionEndReply`, which is an empty message indicating
        successful termination of the interaction.

        :raises ValueError: if the there is no IMD connection available, or
            if invalid parameters are passed to the server.
        :raises KeyError: if the given interaction ID does not exist.

        """
        return self._imd_client.stop_interaction(interaction_id)

    def run_play(self):
        """
        Sends a request to start playing the trajectory to the trajectory service.
        """
        self._frame_client.run_command(PLAY_COMMAND_KEY)

    def run_step(self):
        """
        Sends a request to take one step to the trajectory service.
        """
        self._frame_client.run_command(STEP_COMMAND_KEY)

    def run_pause(self):
        """
        Sends a request to pause the simulation to the trajectory service.
        """
        self._frame_client.run_command(PAUSE_COMMAND_KEY)

    def run_reset(self):
        """
        Sends a request to reset the simulation to the trajectory service.
        """
        self._frame_client.run_command(RESET_COMMAND_KEY)

    def update_available_commands(self) -> MutableMapping[str, CommandInfo]:
        """
        Fetches an updated set of available commands from the services this client is connected
        to.

        :return: A collection of :class:`CommandInfo`, detailing the commands available.

        If the same command name is available on multiple services, the nested nature of the
        returned :class:`ChainMap` will enable the user to determine the correct one to call.
        """

        self._trajectory_commands = _update_commands(self._frame_client)
        self._imd_commands = _update_commands(self._imd_client)
        self._multiplayer_commands = _update_commands(self._multiplayer_client)
        return ChainMap(self._trajectory_commands, self._imd_commands, self._multiplayer_commands)

    def run_command(self, name, **args):
        """
        Runs a command on the trajectory service, multiplayer service or imd service as appropriate.

        :param name: Name of the command to run
        :param args: Dictionary of arguments to run with the command.
        :return: Results of the command, if any.
        """

        if name in self._trajectory_commands:
            return self.run_trajectory_command(name, **args)
        if name in self._imd_commands:
            return self.run_imd_command(name, **args)
        if name in self._multiplayer_commands:
            return self.run_multiplayer_command(name, **args)
        else:
            raise KeyError(f"Unknown command: {name}, run update_available_commands to refresh commands.")

    def run_trajectory_command(self, name: str, **args) -> Dict[str, object]:
        """
        Runs a command on the trajectory service.

        :param name: Name of the command to run
        :param args: Dictionary of arguments to run with the command.
        :return: Results of the command, if any.
        """
        return self._frame_client.run_command(name, **args)

    def run_imd_command(self, name: str, **args) -> Dict[str, object]:
        """
        Runs a command on the iMD service.

        :param name: Name of the command to run
        :param args: Dictionary of arguments to run with the command.
        :return: Results of the command, if any.
        """

        return self._imd_client.run_command(name, **args)

    def run_multiplayer_command(self, name: str, **args):
        """
        Runs a command on the multiplayer service.

        :param name: Name of the command to run
        :param args: Dictionary of arguments to run with the command.
        :return: Results of the command, if any.
        """

        return self._multiplayer_client.run_command(name, **args)

    def join_multiplayer(self, player_name):
        """
        Joins multiplayer with the given player name.

        :param player_name: The player name with which to be identified.

        :raises grpc._channel._Rendezvous: When not connected to a
            multiplayer service
        """
        self._multiplayer_client.join_multiplayer(player_name)

    def set_shared_value(self, key, value) -> bool:
        """
        Attempts to set the given key/value pair on the multiplayer shared value store.

        :param key: The key that identifies the value to be stored.
        :param value: The new value to store.
        :return: `True` if successful, `False` otherwise.


        :raises grpc._channel._Rendezvous: When not connected to a
            multiplayer service
        """
        return self._multiplayer_client.try_set_resource_value(key, value)

    def remove_shared_value(self, key: str) -> bool:
        """
        Attempts to remove the given key on the multiplayer shared value store.

        :raises grpc._channel._Rendezvous: When not connected to a
            multiplayer service
        """
        return self._multiplayer_client.try_remove_resource_key(key)

    def get_shared_value(self, key):
        """
        Attempts to retrieve the value for the given key in the multiplayer shared value store.

        :param key: The key that identifies the value
        :return: The value stored in the dictionary

        :raises grpc._channel._Rendezvous: When not connected to a
            multiplayer service
        """
        return self._multiplayer_client.resources[key]

    @property
    def root_selection(self) -> RenderingSelection:
        """
        Get the root selection, creating it if it does not exist yet.

        :return: The selection representing the root selection of the system
        """
        try:
            root_selection = self.get_selection(SELECTION_ROOT_ID)
        except KeyError:
            root_selection = self._create_selection_from_id_and_name(SELECTION_ROOT_ID, SELECTION_ROOT_NAME)
        root_selection.selected_particle_ids = None
        return root_selection

    def create_selection(
            self,
            name: str,
            particle_ids: Optional[Iterable[int]] = None,
    ) -> RenderingSelection:
        """
        Create a particle selection with the given name.

        :param name: The user-friendly name of the selection.
        :param particle_ids: The indices of the particles to include in the selection.
        :return: The selection that was created.
        """
        if particle_ids is None:
            particle_ids = set()

        # Give the selection an ID based upon the multiplayer player ID and an incrementing counter
        selection_id = f'selection.{self._multiplayer_client.player_id}.{self._next_selection_id}'
        self._next_selection_id += 1

        # Create the selection and setup the particles that it contains
        selection = self._create_selection_from_id_and_name(selection_id, name)
        selection.set_particles(particle_ids)

        # Mark the selection as needing updating, which adds it to the shared value store.
        self.update_selection(selection)

        return selection

    def update_selection(self, selection: RenderingSelection):
        """
        Applies changes to the given selection to the shared key store.

        :param selection: The selection to update.
        """
        struct = dict_to_struct(selection.to_dictionary())
        self.set_shared_value(selection.selection_id, Value(struct_value=struct))

    def remove_selection(self, selection: RenderingSelection):
        """
        Delete the given selection
        """
        self.remove_shared_value(selection.selection_id)

    def clear_selections(self):
        """
        Remove all selections in the system
        """
        selections = list(self.selections)
        for selection in selections:
            self.remove_selection(selection)

    @property
    def selections(self) -> Iterable[RenderingSelection]:
        """
        Get all selections which are stored in the shared key store.

        :return: An iterable of all the selections stored in the shared key store.
        """
        for key, value in self._multiplayer_client.resources.items():
            if key.startswith('selection.'):
                yield self.get_selection(key)

    def get_selection(self, id: str) -> RenderingSelection:
        """
        Get the selection with the given selection id, throwing a KeyError if
        it is not present. For the root selection, use the root_selection
        property.

        :param id: The id of the selection
        :return: The selection if it is present
        """
        value = self._multiplayer_client.resources[id]
        return self._create_selection_from_protobuf_value(value)

    def _create_selection_from_protobuf_value(self, value: Value) -> RenderingSelection:

        selection = RenderingSelection.from_dictionary(struct_to_dict(value.struct_value))
        selection.updated.add_callback(self.update_selection)
        selection.removed.add_callback(self.remove_selection)
        return selection

    def _create_selection_from_id_and_name(self, id: str, name: str) -> RenderingSelection:
        selection = RenderingSelection(id, name)
        selection.updated.add_callback(self.update_selection)
        selection.removed.add_callback(self.remove_selection)
        return selection

    def _join_trajectory(self):
        if self.all_frames:
            self._frame_client.subscribe_frames_async(self._on_frame_received)
        else:
            self._frame_client.subscribe_last_frames_async(
                self._on_frame_received,
                DEFAULT_SUBSCRIPTION_INTERVAL,
            )

    def _join_interactions(self):
        self._imd_client.subscribe_interactions(
            interval=DEFAULT_SUBSCRIPTION_INTERVAL)

    def _on_frame_received(self, frame_index: int, frame: FrameData):
        if self._first_frame is None:
            self._first_frame = frame
        self._frames.append(frame)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
