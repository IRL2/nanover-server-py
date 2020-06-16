# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.
import logging
from concurrent.futures import Future
from queue import Queue
from typing import Iterable, Optional, Dict, Any, NamedTuple, Union

import grpc
from narupa.core import get_requested_port_or_default, NarupaStubClient
from narupa.imd.imd_server import DEFAULT_PORT
from narupa.imd.particle_interaction import ParticleInteraction
from narupa.protocol.imd import (
    InteractiveMolecularDynamicsStub,
    InteractionEndReply,
    SubscribeInteractionsRequest,
)


class Sentinel:
    """
    A specific message passed as a sentinel.
    """


class _LocalInteractionState(NamedTuple):
    """
    Internal state for managing and publishing the client's active interactions.
    """
    queue: Queue
    sentinel: Sentinel
    future: Future


def queue_generator(queue: Queue, sentinel: object):
    """
    Produces a generator that can be used to iterate over the values submitted
    to a queue, until the given sentinel is added to the queue. The iterator
    will block until this sentinel is passed.

    Enables one to take control of when items are passed to a streaming
    iterator, such as that used in :func:`publish_interactions_async`.

    :param queue: Queue that items will be submitted to.
    :param sentinel: A sentinel that indicates the end of iteration. When added
        to the queue, the generator stops.
    :return: Yields the items in put into the queue, in the order they were put in.

    >>> queue = Queue()
    >>> sentinel = Sentinel()
    >>> queue.put(1)
    >>> queue.put(2)
    >>> queue.put(sentinel)
    >>> generator = queue_generator(queue, sentinel)
    >>> for item in generator:
    ...     print(item)
    1
    2

    """
    for val in iter(queue.get, sentinel):
        yield val


def _to_proto(interactions: Iterable[ParticleInteraction]):
    for interaction in interactions:
        yield interaction.proto


class ImdClient(NarupaStubClient):
    """
    A simple IMD client, primarily for testing the IMD server.

    """
    stub: InteractiveMolecularDynamicsStub
    interactions: Dict[str, ParticleInteraction]
    _next_local_interaction_id: int = 0
    _local_interactions_states: Dict[str, _LocalInteractionState]
    _logger: logging.Logger

    def __init__(self, *,
                 channel: grpc.Channel,
                 make_channel_owner: bool = False):
        super().__init__(channel=channel, stub=InteractiveMolecularDynamicsStub, make_channel_owner=make_channel_owner)
        self.interactions = {}
        self._local_interactions_states = {}
        self._logger = logging.getLogger(__name__)

    def start_interaction(self) -> str:
        """
        Start an interaction

        :return: A unique identifier to be used to update the interaction.
        """
        queue: Queue[Union[ParticleInteraction, Sentinel]] = Queue()
        sentinel = Sentinel()
        future = self.publish_interactions_async(queue_generator(queue, sentinel))
        local_interaction_id = self._get_new_local_interaction_id()
        self._local_interactions_states[local_interaction_id] = _LocalInteractionState(queue, sentinel, future)
        return local_interaction_id

    def update_interaction(
            self,
            local_interaction_id: str,
            interaction: ParticleInteraction,
    ):
        """
        Updates the interaction identified with the given interaction_id on the server with
        parameters from the given interaction.

        :param local_interaction_id: The unique interaction ID, created with
            :func:`~ImdClient.start_interaction`, that identifies the
            interaction to update.
        :param interaction: The :class: ParticleInteraction providing new
            parameters for the interaction.

        :raises: ValueError, if invalid parameters are passed to the server.
        :raises: KeyError, if the given interaction ID does not exist.
        """
        if local_interaction_id not in self._local_interactions_states:
            raise KeyError("Attempted to update an interaction with an unknown interaction ID.")
        queue, _, _ = self._local_interactions_states[local_interaction_id]
        queue.put(interaction)

    def stop_interaction(self, local_interaction_id: str) -> InteractionEndReply:
        """
        Stops the interaction identified with the given interaction_id on the server.

        :param local_interaction_id: The unique interaction ID, created with
            :func:`~ImdClient.start_interaction`, that identifies the
            interaction to stop.

        :raises: KeyError, if the given interaction ID does not exist.
        """
        if local_interaction_id not in self._local_interactions_states:
            raise KeyError("Attempted to stop an interaction with an unknown interaction ID.")
        queue, sentinel, future = self._local_interactions_states[local_interaction_id]
        queue.put(sentinel)
        del self._local_interactions_states[local_interaction_id]
        return future.result()

    def publish_interactions_async(self, interactions: Iterable[ParticleInteraction]) -> Future:
        """
        Publishes the iterable of interactions on a thread.

        :param interactions: An iterable generator or collection of interactions.
        :return Future representing the state of the interaction task.
        """
        return self.threads.submit(self.publish_interactions, interactions)

    def publish_interactions(self, interactions: Iterable[ParticleInteraction]) -> InteractionEndReply:
        """
        Publishes the generator of interactions on a thread.

        :param interactions: An iterable generator or collection of interactions.
        :return: A reply indicating successful publishing of interaction.
        """
        return self.stub.PublishInteraction(_to_proto(interactions))

    def stop_all_interactions(self):
        """
        Stops all active interactions governed by this client.
        """
        for local_interaction_id in list(self._local_interactions_states.keys()):
            self.stop_interaction(local_interaction_id)

    def subscribe_interactions(self, interval=0) -> Future:
        """
        Begin receiving updates known interactions.

        :param interval: Minimum time (in seconds) between receiving new updates
            for interaction changes.
        """
        request = SubscribeInteractionsRequest(update_interval=interval)
        return self.threads.submit(self._subscribe_interactions, request)

    def _get_new_local_interaction_id(self) -> str:
        local_interaction_id = str(self._next_local_interaction_id)
        self._next_local_interaction_id += 1
        return local_interaction_id

    def _subscribe_interactions(self, request: SubscribeInteractionsRequest):
        for update in self.stub.SubscribeInteractions(request):
            for interaction_id in update.removals:
                removed = self.interactions.pop(interaction_id, None)
            for interaction in update.updated_interactions:
                self.interactions[interaction.interaction_id] = ParticleInteraction.from_proto(interaction)

    def close(self):
        """
        Closes the IMD client.
        """
        try:
            self.stop_all_interactions()
        except grpc.RpcError as e:
            self._logger.exception(e)
            raise e
        finally:
            super().close()
