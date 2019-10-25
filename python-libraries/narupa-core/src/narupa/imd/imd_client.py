# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.
import logging
import time
from concurrent import futures
from concurrent.futures import Future
from queue import Queue
from typing import Iterable, Optional

import grpc
from narupa.core import get_requested_port_or_default, DEFAULT_CONNECT_ADDRESS, \
    GrpcClient
from narupa.imd.imd_server import DEFAULT_PORT
from narupa.imd.particle_interaction import ParticleInteraction
from narupa.protocol.imd import InteractiveMolecularDynamicsStub, InteractionEndReply


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
    >>> sentinel = object()
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


class ImdClient(GrpcClient):
    """
    A simple IMD client, primarily for testing the IMD server.

    :param address: Address of the IMD server to connect to.
    :param port: The port of the IMD server to connect to.
    """
    def __init__(self, *, address: Optional[str] = None,
                 port: Optional[int] = None):
        port = get_requested_port_or_default(port, DEFAULT_PORT)
        super().__init__(address=address, port=port,
                         stub=InteractiveMolecularDynamicsStub)
        self._active_interactions = {}
        self._logger = logging.getLogger(__name__)

    def start_interaction(self) -> int:
        """
        Start an interaction

        :return: A unique identifier to be used to update the interaction.
        """
        queue = Queue()
        sentinel = object()
        future = self.publish_interactions_async(queue_generator(queue, sentinel))
        if len(self._active_interactions) == 0:
            interaction_id = 0
        else:
            interaction_id = max(self._active_interactions) + 1
        self._active_interactions[interaction_id] = (queue, sentinel, future)
        return interaction_id

    def update_interaction(self, interaction_id, interaction: ParticleInteraction):
        """
        Updates the interaction identified with the given interaction_id on the server with
        parameters from the given interaction.

        :param interaction_id: The unique interaction ID, created with
            :func:`~ImdClient.start_interaction`, that identifies the
            interaction to update.
        :param interaction: The :class: ParticleInteraction providing new
            parameters for the interaction.

        :raises: ValueError, if invalid parameters are passed to the server.
        :raises: KeyError, if the given interaction ID does not exist.
        """
        if interaction_id not in self._active_interactions:
            raise KeyError("Attempted to update an interaction with an unknown interaction ID.")
        queue, _, _ = self._active_interactions[interaction_id]
        queue.put(interaction)

    def stop_interaction(self, interaction_id) -> InteractionEndReply:
        """
        Stops the interaction identified with the given interaction_id on the server.

        :param interaction_id: The unique interaction ID, created with
            :func:`~ImdClient.start_interaction`, that identifies the
            interaction to stop.

        :raises: KeyError, if the given interaction ID does not exist.
        """
        if interaction_id not in self._active_interactions:
            raise KeyError("Attempted to stop an interaction with an unknown interaction ID.")
        queue, sentinel, future = self._active_interactions[interaction_id]
        queue.put(sentinel)
        del self._active_interactions[interaction_id]
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
        for interaction_id in list(self._active_interactions.keys()):
            self.stop_interaction(interaction_id)

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
