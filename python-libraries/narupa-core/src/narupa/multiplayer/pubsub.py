# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.

"""
Module providing a PubSub (publication/subscribe) model for use with gRPC.
"""
import logging
from queue import Queue, Empty
from typing import Dict, Iterator

from narupa.core.request_queues import SingleItemQueue, DictOfQueues

import grpc


class PubSub:
    """
    Represents a multiplayer PubSub service, which simply forwards all publications to all other clients.

    This is achieved with two streams, one for publishing and one for subscribing, and queues.
    The subscription and publication requests must contain a hashable player_id field.

    :param max_queue_size: The maximum number of entries a queue will allow before blocking.
    :param send_latest: Whether to send the latest publication to a client when it subscribes.
    """

    def __init__(self, max_queue_size=10, send_latest=True):
        self.send_latest = send_latest
        self.queues = DictOfQueues()
        self.max_queue_size = max_queue_size
        self.latest_publication = None
        self.logger = logging.Logger(__name__)

    def subscribe(self, request, context):
        """
        Subscribe to a (potentially infinite) stream.

        :param request: Request with player_id field.
        :param context:
        :return: Publications from queue.
        """
        queue = SingleItemQueue(maxsize=self.max_queue_size)
        if request.player_id in self.queues:
            raise KeyError("Player ID already subscribed! Join multiplayer to get a unique player ID.")
        
        with self.queues.one_queue(request.player_id) as queue:
            if self.send_latest and self.latest_publication is not None:
                self.logger.debug(f'Sending latest avatar to: {request.player_id}.')
                yield self.latest_publication
            while context.is_active():
                # wait for the next publication and yield it.
                try:
                    publication = queue.get(block=True, timeout=0.5)
                except Empty:
                    pass
                else:
                    self.latest_publication = publication
                    self.logger.debug(f'Sending avatar to: {request.player_id}.')
                    yield publication

    def start_publish(self, request_iterator: Iterator, context,
                      send_self=False):
        """
        Publish to a stream from a (potentially infinite) iterator.

        :param request_iterator: Request with player_id field.
        :param context:
        :param send_self: Whether to send publications to self (use for debugging).
        :return:
        """
        try:
            # TODO no clean way to elegantly cancel publication?
            for request in request_iterator:
                self.publish(request, send_self)
        except grpc.RpcError as e:
            if not context.is_active():
                return
            else:  # pragma: no cover
                raise e  # pragma: no cover

    def publish(self, publication, send_self=False):
        """
        Publish an object onto the stream.

        :param publication: A new publication with player_id field.
        :param send_self: Whether to send publications to self (use for debugging).
        :return:
        """
        for player_id, queue in self.queues.iter_queues_items():
            if not send_self and (publication.player_id == player_id):
                continue
            queue.put(publication)
