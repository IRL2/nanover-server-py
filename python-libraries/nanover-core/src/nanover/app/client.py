"""
Module containing a basic interactive molecular dynamics client that receives frames
and can publish interactions.
"""

import time
import warnings
from collections import deque, ChainMap
from functools import wraps, partial
from typing import Iterable, Tuple, Type, TypeVar, cast
from typing import Optional, Sequence, Dict, MutableMapping
from uuid import uuid4

import grpc
from grpc import RpcError, StatusCode
from nanover.app.app_server import DEFAULT_NANOVER_PORT, MULTIPLAYER_SERVICE_NAME
from nanover.app.selection import RenderingSelection
from nanover.command import CommandInfo
from nanover.core import NanoverClient, DEFAULT_CONNECT_ADDRESS
from nanover.core.nanover_client import DEFAULT_STATE_UPDATE_INTERVAL
from nanover.essd import DiscoveryClient
from nanover.imd import ImdClient, IMD_SERVICE_NAME
from nanover.imd.particle_interaction import ParticleInteraction
from nanover.trajectory import FrameClient, FrameData, FRAME_SERVICE_NAME
from nanover.trajectory.frame_server import (
    PLAY_COMMAND_KEY,
    STEP_COMMAND_KEY,
    PAUSE_COMMAND_KEY,
    RESET_COMMAND_KEY,
)
from nanover.utilities.change_buffers import DictionaryChange

# Default to a low framerate to avoid build up in the frame stream
DEFAULT_SUBSCRIPTION_INTERVAL = 1 / 30

# ID of the root selection
SELECTION_ROOT_ID = "selection.root"
# Name of the root selection
SELECTION_ROOT_NAME = "Root Selection"


ClientVarType = TypeVar("ClientVarType", bound=NanoverClient)


def _update_commands(client: Optional[NanoverClient]):
    if client is None:
        return {}
    try:
        return client.update_available_commands()
    except RpcError as e:
        if e._state.code == StatusCode.UNAVAILABLE:
            return {}
        raise e


def _need_attribute(func, *, name: str, attr: str):
    @wraps(func)
    def wrapper(self, *args, **kwargs):
        if getattr(self, attr) is None:
            raise RuntimeError(f"Not connected to {name} service")
        return func(self, *args, **kwargs)

    return wrapper


def need_trajectory_joined(func):
    @wraps(func)
    def wrapper(self, *args, **kwargs):
        if not self._are_framed_subscribed:
            raise RuntimeError(
                "You need to first subscribe to the frame using the "
                "subscribe_to_frames or the subscribe_to_all_frames methods."
            )
        return func(self, *args, **kwargs)

    return wrapper


# Use partial to specify which attribute is needed for the given decorator.
need_frames = partial(_need_attribute, name="trajectory", attr="_frame_client")
need_imd = partial(_need_attribute, name="imd", attr="_imd_client")
need_multiplayer = partial(
    _need_attribute, name="multiplayer", attr="_multiplayer_client"
)


class NanoverImdClient:
    """
    Interactive molecular dynamics client that receives frames, create selections,
    and join the multiplayer shared state.

    :param trajectory_address: Address and port of the trajectory service.
    :param imd_address: Address and port of the iMD service.
    :param multiplayer_address: Address and port of the multiplayer service.
    :param max_frames: Maximum number of frames to store in a buffer, if not storing all frames.

    All addresses are optional, so one can, for example, just connect to a trajectory service to passively receive
    frames.

    The :func:`NanoverImdClient.autoconnect` and :func:`NanoverImdClient.connect_to_single_server` methods provide
    shorthands for common server setups.

    Inspecting a Frame
    ==================

    The NanoVer Imd client can be used to inspect frames received from a :class:`nanover.trajectory.FrameServer`,
    which can be useful for analysis.

    .. python
        # Assuming there is only one server (or set of servers) running.
        client = NanoverImdClient.autoconnect()
        # Request to receive the frames from the server
        client.subscribe_to_frames()
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
        client.subscribe_multiplayer()
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

    _player_id: str
    _channels: Dict[Tuple[str, int], grpc.Channel]

    _frame_client: Optional[FrameClient]
    _imd_client: Optional[ImdClient]
    _multiplayer_client: Optional[NanoverClient]
    _frames: deque
    _current_frame: FrameData
    _first_frame: Optional[FrameData]

    _next_selection_id: int = 0

    _trajectory_commands: Dict[str, CommandInfo]
    _imd_commands: Dict[str, CommandInfo]
    _multiplayer_commands: Dict[str, CommandInfo]

    _are_framed_subscribed: bool
    _subscribed_to_all_frames: bool

    def __init__(
        self,
        *,
        trajectory_address: Optional[Tuple[str, int]] = None,
        imd_address: Optional[Tuple[str, int]] = None,
        multiplayer_address: Optional[Tuple[str, int]] = None,
        max_frames=50,
        all_frames: Optional[bool] = None,
    ):
        if all_frames is not None:
            warnings.warn(
                "The `all_frames` argument is deprecated and will be removed "
                "in a later version. Use `subscribe_to_all_frames` to "
                "subscribe to all frames, or `subscribe_to_frames` to "
                "subscribe with a set interval.",
                DeprecationWarning,
            )
        self._player_id = str(uuid4())

        self._channels = {}

        self.max_frames = max_frames

        self._frame_client = None
        self._multiplayer_client = None
        self._imd_client = None
        self.connect(
            trajectory_address=trajectory_address,
            imd_address=imd_address,
            multiplayer_address=multiplayer_address,
        )

        self._are_framed_subscribed = False
        self._subscribed_to_all_frames = False
        self._frames = deque(maxlen=self.max_frames)
        self._first_frame = None
        self._current_frame = FrameData()

        # initialise the set of available commands.
        self.update_available_commands()

    @property
    def all_frames(self):
        return self._subscribed_to_all_frames

    @classmethod
    def connect_to_single_server(
        cls, address: Optional[str] = None, port: Optional[int] = None
    ):
        """
        Connect to a single NanoVer server running all services on the same port.

        :param address: Address of the server.
        :param port: Server port
        :return: Instantiation of a client connected to all available services on the server at the given destination.
        """
        address = address or DEFAULT_CONNECT_ADDRESS
        port = port or DEFAULT_NANOVER_PORT
        url = (address, port)
        return cls(trajectory_address=url, imd_address=url, multiplayer_address=url)

    @classmethod
    def connect_to_single_server_multiple_ports(
        cls,
        address: str,
        trajectory_port: int,
        imd_port: int,
        multiplayer_port: int,
    ):
        """
        Connect to a collection of NanoVer servers running at the same address but potentially different ports.

        :param address: Address of the server.
        :param multiplayer_port: The port at which multiplayer is running.
        :param trajectory_port: The port at which the trajectory service is running.
        :param imd_port: The port at which the iMD service is running.
        :return: Instantiation of a client connected to all available services on the server at the given destination.
        """

        # TODO this is a utility method for testing... a good place to put this?
        return cls(
            trajectory_address=(address, trajectory_port),
            imd_address=(address, imd_port),
            multiplayer_address=(address, multiplayer_port),
        )

    @classmethod
    def autoconnect(
        cls,
        search_time=2.0,
        discovery_address: Optional[str] = None,
        discovery_port: Optional[int] = None,
        name: Optional[str] = None,
    ):
        """
        Autoconnect to the first available server discovered that at least produces frames.

        :param search_time: Time, in seconds, to search for.
        :param discovery_address: IP address to search on.
        :param discovery_port: Port upon which to listen for discovery messages.
        :param name: If supplied, only servers with this name will be used.
        :return: Instantiation of an iMD client connected to whatever is available at the first
        """
        if name is not None:
            first_service = _search_for_first_server_with_name(
                name, search_time, discovery_address, discovery_port
            )
        else:
            first_service = _search_for_first_available_frame_service(
                search_time, discovery_address, discovery_port
            )

        if first_service is None:
            raise ConnectionError("Could not find an iMD server")

        trajectory_address = first_service.get_service_address(FRAME_SERVICE_NAME)
        imd_address = first_service.get_service_address(IMD_SERVICE_NAME)
        multiplayer_address = first_service.get_service_address(
            MULTIPLAYER_SERVICE_NAME
        )
        return cls(
            trajectory_address=trajectory_address,
            imd_address=imd_address,
            multiplayer_address=multiplayer_address,
        )

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
            # We are not subscribed to the frames anymore.
            self._are_framed_subscribed = False
            self._subscribed_to_all_frames = False
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
        self._multiplayer_client = self._connect_client(NanoverClient, address)

    def connect(
        self,
        *,
        trajectory_address: Optional[Tuple[str, int]] = None,
        imd_address: Optional[Tuple[str, int]] = None,
        multiplayer_address: Optional[Tuple[str, int]] = None,
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
    @need_trajectory_joined
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

    @property  # type: ignore
    @need_frames
    @need_trajectory_joined
    def latest_frame(self) -> Optional[FrameData]:
        """
        The trajectory frame most recently received, if any.

        :return: :class:`FrameData`, or `None` if none has been received.
        """
        if len(self.frames) == 0:
            return None
        return self.frames[-1]

    @property  # type: ignore
    @need_frames
    @need_trajectory_joined
    def current_frame(self) -> FrameData:
        """
        A copy of the current state of the trajectory, formed by collating all received frames.

        :return: :class:`FrameData`, which is empty if none have been received.
        """
        return self._current_frame.copy()

    @property  # type: ignore
    @need_multiplayer
    def latest_multiplayer_values(self) -> Dict[str, object]:
        """
        The latest state of the multiplayer shared key/value store.

        :return: Dictionary of the current state of multiplayer shared key/value store.
        """
        return self._multiplayer_client.copy_state()  # type: ignore

    @property  # type: ignore
    @need_frames
    @need_trajectory_joined
    def frames(self) -> Sequence[FrameData]:
        """
        The most recently received frames up to the storage limit specified
        by `max_frames`.

        :return: Sequence of frames.
        """
        return list(self._frames)

    @property  # type: ignore
    @need_frames
    @need_trajectory_joined
    def first_frame(self) -> Optional[FrameData]:
        """
        The first received trajectory frame, if any.

        :return: The first frame received by this trajectory, or `None`.
        """
        return self._first_frame

    @property  # type: ignore
    @need_imd
    def interactions(self) -> Dict[str, ParticleInteraction]:
        """
        The dictionary of current interactions received by this client.
        :return: Dictionary of active interactions, keyed by interaction ID identifying who is performing the
        interactions.
        """
        return self._imd_client.interactions  # type: ignore

    @need_imd
    def start_interaction(
        self, interaction: Optional[ParticleInteraction] = None
    ) -> str:
        """
        Start an interaction with the IMD server.

        :param interaction: An optional :class: ParticleInteraction with which
            to begin.
        :return: The unique interaction ID of this interaction, which can be
            used to update the interaction with
            :func:`~NanoverClient.update_interaction`.

        :raises: ValueError, if the there is no IMD connection available.
        """
        interaction_id = self._imd_client.start_interaction()  # type: ignore
        if interaction is not None:
            self.update_interaction(interaction_id, interaction)
        return interaction_id

    @need_imd
    def update_interaction(self, interaction_id, interaction: ParticleInteraction):
        """
        Updates the interaction identified with the given interaction_id on
        the server with parameters from the given interaction.

        :param interaction_id: The unique id of the interaction to be updated.
        :param interaction: The :class: ParticleInteraction providing new
            parameters for the interaction.

        :raises: ValueError, if the there is no IMD connection available, or
            if invalid parameters are passed to the server.
        :raises KeyError: if the given interaction ID does not exist.
        """
        self._imd_client.update_interaction(interaction_id, interaction)  # type: ignore

    @need_imd
    def stop_interaction(self, interaction_id) -> bool:
        """
        Stops the interaction identified with the given interaction_id on the server.
        This method blocks until the server returns a reply indicating that the
        interaction has stopped.

        :param interaction_id: The unique interaction ID, created with
            :func:`~NanoverClient.start_interaction`, that identifies the
            interaction to stop.
        :return: An :class:`InteractionEndReply`, which is an empty message indicating
            successful termination of the interaction.

        :raises ValueError: if the there is no IMD connection available, or
            if invalid parameters are passed to the server.
        :raises KeyError: if the given interaction ID does not exist.
        """
        return self._imd_client.stop_interaction(interaction_id)  # type: ignore

    @need_frames
    def run_play(self):
        """
        Sends a request to start playing the trajectory to the trajectory service.
        """
        self._frame_client.run_command(PLAY_COMMAND_KEY)  # type: ignore

    @need_frames
    def run_step(self):
        """
        Sends a request to take one step to the trajectory service.
        """
        self._frame_client.run_command(STEP_COMMAND_KEY)  # type: ignore

    @need_frames
    def run_pause(self):
        """
        Sends a request to pause the simulation to the trajectory service.
        """
        self._frame_client.run_command(PAUSE_COMMAND_KEY)  # type: ignore

    @need_frames
    def run_reset(self):
        """
        Sends a request to reset the simulation to the trajectory service.
        """
        self._frame_client.run_command(RESET_COMMAND_KEY)  # type: ignore

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
        return ChainMap(
            self._trajectory_commands, self._imd_commands, self._multiplayer_commands
        )

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
        raise KeyError(
            f"Unknown command: {name}, run update_available_commands to refresh commands."
        )

    @need_frames
    def run_trajectory_command(self, name: str, **args) -> Dict[str, object]:
        """
        Runs a command on the trajectory service.

        :param name: Name of the command to run
        :param args: Dictionary of arguments to run with the command.
        :return: Results of the command, if any.
        """
        return self._frame_client.run_command(name, **args)  # type: ignore

    @need_imd
    def run_imd_command(self, name: str, **args) -> Dict[str, object]:
        """
        Runs a command on the iMD service.

        :param name: Name of the command to run
        :param args: Dictionary of arguments to run with the command.
        :return: Results of the command, if any.
        """

        return self._imd_client.run_command(name, **args)  # type: ignore

    @need_multiplayer
    def run_multiplayer_command(self, name: str, **args):
        """
        Runs a command on the multiplayer service.

        :param name: Name of the command to run.
        :param args: Dictionary of arguments to run with the command.
        :return: Results of the command, if any.
        """

        return self._multiplayer_client.run_command(name, **args)  # type: ignore

    @need_multiplayer
    def subscribe_multiplayer(self, interval=DEFAULT_STATE_UPDATE_INTERVAL):
        """
        Subscribe to all multiplayer state updates.

        :param interval: Subscription interval for state updates.

        :raises grpc._channel._Rendezvous: When not connected to a
            multiplayer service
        """
        multiplayer_client = cast(NanoverClient, self._multiplayer_client)
        multiplayer_client.subscribe_all_state_updates(interval)

    @need_multiplayer
    def attempt_update_multiplayer_state(
        self,
        update: DictionaryChange,
    ) -> bool:
        """
        Attempt to make a single atomic change to the shared state, blocking
        until a response is received.

        :param update: A single change to make to the shared state that will
            either be made in full, or ignored if some of the keys are locked
            by another user.
        :return: True if the server accepted our change, and False otherwise.
        """
        return self._multiplayer_client.attempt_update_state(update)  # type: ignore

    @need_multiplayer
    def attempt_update_multiplayer_locks(
        self,
        update: Dict[str, Optional[float]],
    ) -> bool:
        """
        Attempt to acquire and/or free a number of locks on the shared state.

        :param update: A dictionary of keys to either a duration in
            seconds to attempt to acquire or renew a lock, or None to indicate
            the lock should be released if held.
        :return: True if the desired locks were acquired, and False otherwise.
        """
        return self._multiplayer_client.attempt_update_locks(update)  # type: ignore

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
        change = DictionaryChange(updates={key: value})
        return self.attempt_update_multiplayer_state(change)  # type: ignore

    @need_multiplayer
    def remove_shared_value(self, key: str) -> bool:
        """
        Attempts to remove the given key on the multiplayer shared value store.

        :raises grpc._channel._Rendezvous: When not connected to a
            multiplayer service
        """
        change = DictionaryChange(removals=set([key]))
        return self.attempt_update_multiplayer_state(change)  # type: ignore

    @need_multiplayer
    def get_shared_value(self, key):
        """
        Attempts to retrieve the value for the given key in the multiplayer shared value store.

        :param key: The key that identifies the value
        :return: The value stored in the dictionary

        :raises grpc._channel._Rendezvous: When not connected to a
            multiplayer service
        """
        return self._multiplayer_client.copy_state()[key]  # type: ignore

    @property  # type: ignore
    @need_multiplayer
    def root_selection(self) -> RenderingSelection:
        """
        Get the root selection, creating it if it does not exist yet.

        :return: The selection representing the root selection of the system
        """
        try:
            root_selection = self.get_selection(SELECTION_ROOT_ID)
        except KeyError:
            root_selection = self._create_selection_from_id_and_name(
                SELECTION_ROOT_ID, SELECTION_ROOT_NAME
            )
        root_selection.selected_particle_ids = set()
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
        selection_id = f"selection.{self._player_id}.{self._next_selection_id}"
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

    @property  # type: ignore
    @need_multiplayer
    def selections(self) -> Iterable[RenderingSelection]:
        """
        Get all selections which are stored in the shared key store.

        :return: An iterable of all the selections stored in the shared key store.
        """
        for key, _ in self._multiplayer_client.copy_state().items():  # type: ignore
            if key.startswith("selection."):
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
        value = self._multiplayer_client.copy_state()[selection_id]  # type: ignore
        return self._create_selection_from_dict(value)

    def _create_selection_from_dict(self, value) -> RenderingSelection:
        selection = RenderingSelection.from_dictionary(value)
        selection.updated.add_callback(self.update_selection)
        selection.removed.add_callback(self.remove_selection)
        return selection

    def _create_selection_from_id_and_name(
        self, selection_id: str, name: str
    ) -> RenderingSelection:
        selection = RenderingSelection(selection_id, name)
        selection.updated.add_callback(self.update_selection)
        selection.removed.add_callback(self.remove_selection)
        return selection

    @property
    def are_frames_subscribed(self):
        """
        Returns `True` if the client is subscribed to frames with either
        :meth:`subscribe_to_all_frames` or :meth:`subscribe_to_frames`.
        """
        return self._are_framed_subscribed

    @need_frames
    def subscribe_to_all_frames(self):
        """
        Request all the frames from the server.

        This makes the frames available with the :attr:`frames`,
        :attr:`current_frame`, :attr:`latest_frames`, and :attr:`first_frame`
        attributes.

        All the frames produced by the frame server are sent from the time of
        the subscription get received by the client. Note that, depending on
        the server set up, this may not be all the frames generated by the
        MD engine. Subscribing to all frames is usefull for some applications
        such as saving a trajectory or analyses, however it may get behind the
        server if the frames are not received and treated as fast as they
        are sent.

        .. warning::

            A client can subscribe to frames only ones and cannot change how
            it subscribes.

        .. seealso:: subscribe_to_frames, are_frames_subscribed

        """
        if self._are_framed_subscribed:
            return
        self._subscribed_to_all_frames = True
        frame_client = cast(
            FrameClient, self._frame_client
        )  # @need_frames makes sure of that
        frame_client.subscribe_frames_async(self._on_frame_received)
        self._are_framed_subscribed = True

    @need_frames
    def subscribe_to_frames(self, interval: float = DEFAULT_SUBSCRIPTION_INTERVAL):
        """
        Request the latest frames at a given time interval.

        This makes the frames available with the :attr:`frames`,
        :attr:`current_frame`, :attr:`latest_frames`, and :attr:`first_frame`
        attributes.

        This requests for the server to send the latest available frame at a
        given framerate expressed with the interval of time between two
        consecutive frames. This requested interval is the a minimum interval:
        if the frames are produced slower, then they will be sent as fast as
        possible.

        This is usefull for applications that do not need all the frames but
        only the latest available ones.

        .. warning::

            A client can subscribe to frames only ones and cannot change how
            it subscribes.

        .. seealso:: subscribe_to_all_frames, are_frames_subscribed

        :param interval: Minimum time, in seconds, between two consecutive
            frames.
        """
        if self._are_framed_subscribed:
            return
        self._subscribed_to_all_frames = False
        frame_client = cast(
            FrameClient, self._frame_client
        )  # @need_frames makes sure of that
        frame_client.subscribe_last_frames_async(
            self._on_frame_received,
            DEFAULT_SUBSCRIPTION_INTERVAL,
        )
        self._are_framed_subscribed = True

    def _on_frame_received(self, frame_index: int, frame: FrameData):
        if self._first_frame is None:
            self._first_frame = frame
        self._frames.append(frame)
        self._current_frame.raw.MergeFrom(frame.raw)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def _connect_client(
        self,
        client_type: Type[ClientVarType],
        address: Tuple[str, int],
    ) -> ClientVarType:
        # TODO add support for encryption here somehow.

        # if there already exists a channel with the same address, reuse it, otherwise create a new insecure
        # connection.
        if address in self._channels:
            client: ClientVarType = client_type(
                channel=self._channels[address], make_channel_owner=False
            )
        else:
            client: ClientVarType = client_type.insecure_channel(address=address[0], port=address[1])  # type: ignore[no-redef]
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
