# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.
"""
Module containing a basic interactive molecular dynamics client that receives frames
and can publish interactions.
"""
import time
from collections import deque, ChainMap
from functools import wraps, partial
from typing import Iterable, Tuple, Type
from typing import Optional, Sequence, Dict, MutableMapping

from grpc import RpcError, StatusCode
from narupa.app.app_server import DEFAULT_NARUPA_PORT
from narupa.app.selection import RenderingSelection
from narupa.command import CommandInfo
from narupa.core import NarupaClient, DEFAULT_CONNECT_ADDRESS
from narupa.essd import DiscoveryClient
from narupa.imd import ImdClient, IMD_SERVICE_NAME
from narupa.imd.particle_interaction import ParticleInteraction
from narupa.multiplayer import MultiplayerClient, MULTIPLAYER_SERVICE_NAME
from narupa.protocol.imd import InteractionEndReply
from narupa.trajectory import FrameClient, FrameData, FRAME_SERVICE_NAME
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
    except AttributeError:
        return {}
    except RpcError as e:
        if e._state.code == StatusCode.UNAVAILABLE:
            return {}
        raise e


def _need_attribute(func, *, name: str, attr: str):
    @wraps(func)
    def wrapper(self, *args, **kwargs):
        if getattr(self, attr) is None:
            raise RuntimeError(f'Not connected to {name} service')
        return func(self, *args, **kwargs)

    return wrapper


# Use partial to specify which attribute is needed for the given decorator.
need_frames = partial(_need_attribute, name='trajectory', attr='_frame_client')
need_imd = partial(_need_attribute, name='imd', attr='_imd_client')
need_multiplayer = partial(_need_attribute, name='multiplayer', attr='_multiplayer_client')


class NarupaImdClient:
    """
    Interactive molecular dynamics client that receives frames, create selections,
    and join the multiplayer shared state.

    :param trajectory_address: Address and port of the trajectory service.
    :param imd_address: Address and port of the iMD service.
    :param multiplayer_address: Address and port of the multiplayer service.
    :param max_frames: Maximum number of frames to store in a buffer, if not storing all frames.
    :param all_frames: Whether to receive all frames, or skip to the latest one.

    All addresses are optional, so one can, for example, just connect to a trajectory service to passively receive
    frames.

    The :fun:`NarupaImdClient.autoconnect` and :fun:`NarupaImdClient.connect_to_single_server` methods provide
    shorthands for common server setups.

    Inspecting a Frame
    ==================

    The Narupa Imd client can be used to inspect frames received from a :class:`narupa.trajectory.FrameServer`,
    which can be useful for analysis.

    .. python
        # Assuming there is only one server (or set of servers) running.
        client = NarupaImdClient.autoconnect()
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
    _frame_client: Optional[FrameClient]
    _imd_client: Optional[ImdClient]
    _multiplayer_client: Optional[MultiplayerClient]
    _frames: deque

    _next_selection_id: int = 0

    def __init__(self, *,
                 trajectory_address: Tuple[str, int] = None,
                 imd_address: Tuple[str, int] = None,
                 multiplayer_address: Tuple[str, int] = None,
                 max_frames=50,
                 all_frames=True):

        self._channels = {}

        self.all_frames = all_frames
        self.max_frames = max_frames

        self._frame_client = None
        self._multiplayer_client = None
        self._imd_client = None
        self.connect(trajectory_address=trajectory_address,
                     imd_address=imd_address,
                     multiplayer_address=multiplayer_address)

        self._frames = deque(maxlen=self.max_frames)
        self._first_frame = None

        self.update_available_commands()  # initialise the set of available commands.

    @classmethod
    def connect_to_single_server(cls, address: Optional[str] = None, port: Optional[int] = None):
        """
        Connect to a single Narupa server running all services on the same port.

        :param address: Address of the server.
        :param port: Server port
        :return: Instantiation of a client connected to all available services on the server at the given destination.
        """
        address = address or DEFAULT_CONNECT_ADDRESS
        port = port or DEFAULT_NARUPA_PORT
        url = (address, port)
        return cls(trajectory_address=url, imd_address=url, multiplayer_address=url)

    @classmethod
    def connect_to_single_server_multiple_ports(cls,
                                                address: Optional[str] = None,
                                                trajectory_port: Optional[int] = None,
                                                imd_port: Optional[int] = None,
                                                multiplayer_port: Optional[int] = None,
                                                ):
        """
        Connect to a collection of Narupa servers running at the same address but potentially different ports.

        :param address: Address of the server.
        :param multiplayer_port: The port at which multiplayer is running.
        :param trajectory_port: The port at which the trajectory service is running.
        :param imd_port: The port at which the iMD service is running.
        :return: Instantiation of a client connected to all available services on the server at the given destination.
        """

        # TODO this is a utility method for testing... a good place to put this?
        address = address or DEFAULT_CONNECT_ADDRESS
        return cls(trajectory_address=(address, trajectory_port),
                   imd_address=(address, imd_port),
                   multiplayer_address=(address, multiplayer_port))

    @classmethod
    def autoconnect(cls, search_time=2.0,
                    discovery_address: Optional[str] = None,
                    discovery_port: Optional[int] = None,
                    name: Optional[str] = None):
        """
        Autoconnect to the first available server discovered that at least produces frames.
        
        :param search_time: Time, in seconds, to search for.
        :param discovery_address: IP address to search on.
        :param discovery_port: Port upon which to listen for discovery messages.
        :param name: If supplied, only servers with this name will be used.
        :return: Instantiation of an iMD client connected to whatever is available at the first
        """
        if name is not None:
            first_service = _search_for_first_server_with_name(name, search_time, discovery_address, discovery_port)
        else:
            first_service = _search_for_first_available_frame_service(search_time, discovery_address, discovery_port)

        if first_service is None:
            raise ConnectionError("Could not find an iMD server")

        trajectory_address = first_service.get_service_address(FRAME_SERVICE_NAME)
        imd_address = first_service.get_service_address(IMD_SERVICE_NAME)
        multiplayer_address = first_service.get_service_address(MULTIPLAYER_SERVICE_NAME)
        return cls(trajectory_address=trajectory_address, imd_address=imd_address,
                   multiplayer_address=multiplayer_address)

    def close(self, clear_frames=True):
        """
        Closes the connection with the server.

        :param clear_frames: Whether to clear the frames received by the
            client, or keep them.
        """
        if self._imd_client is not None:
            self._imd_client.close()
            self._imd_client = None
        if self._multiplayer_client is not None:
            self._multiplayer_client.close()
            self._multiplayer_client = None
        if self._frame_client is not None:
            self._frame_client.close()
            self._frame_client = None
        self._channels.clear()

        if clear_frames:
            self._first_frame = None
            self._frames.clear()

    def connect_trajectory(self, address: Tuple[str, int]):
        """
        Connects the client to the given trajectory server, and begin receiving frames.

        :param address: The address and port of the trajectory server.
        """

        self._frame_client = self._connect_client(FrameClient, address)
        self._join_trajectory()

    def connect_imd(self, address: Tuple[str, int]):
        """
        Connects the client to the given interactive molecular dynamics server,
        allowing it to start publishing interactions.

        :param address: The address and port of the IMD server.
        :param port: The port of the IMD server.
        """
        self._imd_client = self._connect_client(ImdClient, address)
        self._imd_client.subscribe_all_state_updates()
    
    def connect_multiplayer(self, address: Tuple[str, int]):
        """
        Connects the client to the given multiplayer server.

        :param address: The address and port of the multiplayer server.
        """
        self._multiplayer_client = self._connect_client(MultiplayerClient, address)

    def connect(self, *,
                trajectory_address: Tuple[str, int] = None,
                imd_address: Tuple[str, int] = None,
                multiplayer_address: Tuple[str, int] = None,
                ):
        """
        Connects the client to all services for which addresses are provided.
        """
        if trajectory_address is not None:
            self.connect_trajectory(trajectory_address)
        if imd_address is not None:
            self.connect_imd(imd_address)
        if multiplayer_address is not None:
            self.connect_multiplayer(multiplayer_address)

    @need_frames
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
    @need_frames
    def latest_frame(self) -> Optional[FrameData]:
        """
        The trajectory frame most recently received, if any.

        :return: :class:`FrameData`, or `None` if none has been received.
        """
        if len(self.frames) is 0:
            return None
        return self.frames[-1]

    @property
    @need_multiplayer
    def latest_multiplayer_values(self) -> Dict[str, object]:
        """
        The latest state of the multiplayer shared key/value store.

        :return: Dictionary of the current state of multiplayer shared key/value store.
        """
        return dict(self._multiplayer_client.resources)

    @property
    @need_frames
    def frames(self) -> Sequence[FrameData]:
        """
        The most recently received frames up to the storage limit specified
        by `max_frames`.

        :return: Sequence of frames.
        """
        return list(self._frames)

    @property
    @need_frames
    def first_frame(self) -> Optional[FrameData]:
        """
        The first received trajectory frame, if any.

        :return: The first frame received by this trajectory, or `None`.
        """
        return self._first_frame

    @property
    @need_imd
    def interactions(self) -> Dict[str, ParticleInteraction]:
        """
        The dictionary of current interactions received by this client.
        :return: Dictionary of active interactions, keyed by interaction ID identifying who is performing the
        interactions.
        """
        return self._imd_client.interactions

    @need_imd
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

    @need_imd
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

    @need_imd
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

    @need_frames
    def run_play(self):
        """
        Sends a request to start playing the trajectory to the trajectory service.
        """
        self._frame_client.run_command(PLAY_COMMAND_KEY)

    @need_frames
    def run_step(self):
        """
        Sends a request to take one step to the trajectory service.
        """
        self._frame_client.run_command(STEP_COMMAND_KEY)

    @need_frames
    def run_pause(self):
        """
        Sends a request to pause the simulation to the trajectory service.
        """
        self._frame_client.run_command(PAUSE_COMMAND_KEY)

    @need_frames
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
        raise KeyError(f"Unknown command: {name}, run update_available_commands to refresh commands.")

    @need_frames
    def run_trajectory_command(self, name: str, **args) -> Dict[str, object]:
        """
        Runs a command on the trajectory service.

        :param name: Name of the command to run
        :param args: Dictionary of arguments to run with the command.
        :return: Results of the command, if any.
        """
        return self._frame_client.run_command(name, **args)

    @need_imd
    def run_imd_command(self, name: str, **args) -> Dict[str, object]:
        """
        Runs a command on the iMD service.

        :param name: Name of the command to run
        :param args: Dictionary of arguments to run with the command.
        :return: Results of the command, if any.
        """

        return self._imd_client.run_command(name, **args)

    @need_multiplayer
    def run_multiplayer_command(self, name: str, **args):
        """
        Runs a command on the multiplayer service.

        :param name: Name of the command to run
        :param args: Dictionary of arguments to run with the command.
        :return: Results of the command, if any.
        """

        return self._multiplayer_client.run_command(name, **args)

    @need_multiplayer
    def join_multiplayer(self, player_name):
        """
        Joins multiplayer with the given player name.

        :param player_name: The player name with which to be identified.

        :raises grpc._channel._Rendezvous: When not connected to a
            multiplayer service
        """
        self._multiplayer_client.join_multiplayer(player_name)

    @need_multiplayer
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

    @need_multiplayer
    def remove_shared_value(self, key: str) -> bool:
        """
        Attempts to remove the given key on the multiplayer shared value store.

        :raises grpc._channel._Rendezvous: When not connected to a
            multiplayer service
        """
        return self._multiplayer_client.try_remove_resource_key(key)

    @need_multiplayer
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
    @need_multiplayer
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

    @need_multiplayer
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

    @need_multiplayer
    def update_selection(self, selection: RenderingSelection):
        """
        Applies changes to the given selection to the shared key store.

        :param selection: The selection to update.
        """
        self.set_shared_value(selection.selection_id, selection.to_dictionary())

    @need_multiplayer
    def remove_selection(self, selection: RenderingSelection):
        """
        Delete the given selection
        """
        self.remove_shared_value(selection.selection_id)

    @need_multiplayer
    def clear_selections(self):
        """
        Remove all selections in the system
        """
        selections = list(self.selections)
        for selection in selections:
            self.remove_selection(selection)

    @property
    @need_multiplayer
    def selections(self) -> Iterable[RenderingSelection]:
        """
        Get all selections which are stored in the shared key store.

        :return: An iterable of all the selections stored in the shared key store.
        """
        for key, _ in self._multiplayer_client.resources.items():
            if key.startswith('selection.'):
                yield self.get_selection(key)

    @need_multiplayer
    def get_selection(self, selection_id: str) -> RenderingSelection:
        """
        Get the selection with the given selection id, throwing a KeyError if
        it is not present. For the root selection, use the root_selection
        property.

        :param selection_id: The id of the selection
        :return: The selection if it is present
        """
        value = self._multiplayer_client.resources[selection_id]
        return self._create_selection_from_dict(value)

    def _create_selection_from_dict(self, value) -> RenderingSelection:

        selection = RenderingSelection.from_dictionary(value)
        selection.updated.add_callback(self.update_selection)
        selection.removed.add_callback(self.remove_selection)
        return selection

    def _create_selection_from_id_and_name(self, selection_id: str, name: str) -> RenderingSelection:
        selection = RenderingSelection(selection_id, name)
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

    def _on_frame_received(self, frame_index: int, frame: FrameData):
        if self._first_frame is None:
            self._first_frame = frame
        self._frames.append(frame)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def _connect_client(self, client_type: Type[NarupaClient], address: Tuple[str, int]):
        # TODO add support for encryption here somehow.

        # if there already exists a channel with the same address, reuse it, otherwise create a new insecure
        # connection.
        if address in self._channels:
            client = client_type(channel=self._channels[address], make_channel_owner=False)
        else:
            client = client_type.insecure_channel(address=address[0], port=address[1])
            self._channels[address] = client.channel
        return client


def _search_for_first_server_with_name(
        server_name: str,
        search_time: float = 2.0,
        discovery_address: Optional[str] = None,
        discovery_port: Optional[int] = None,
):
    with DiscoveryClient(discovery_address, discovery_port) as discovery_client:
        for hub in discovery_client.search_for_services(search_time):
            if hub.name == server_name:
                return hub
    return None


def _search_for_first_available_frame_service(
        search_time: float = 2.0,
        discovery_address: Optional[str] = None,
        discovery_port: Optional[int] = None,
):
    with DiscoveryClient(discovery_address, discovery_port) as discovery_client:
        for hub in discovery_client.search_for_services(search_time):
            if FRAME_SERVICE_NAME in hub.services:
                return hub
    return None
