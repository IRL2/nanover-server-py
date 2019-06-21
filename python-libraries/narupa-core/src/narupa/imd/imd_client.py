# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.

import time
from concurrent import futures
from concurrent.futures import Future
from typing import Collection, Iterable, Generator

import grpc

from narupa.protocol.imd import InteractiveMolecularDynamicsStub, InteractionEndReply


def delayed_generator(iterable: Iterable, delay: float = 0):
    """
    Turns an iterable collection into a generator, where each item is yielded after the given delay (seconds).
    :param iterable: Iterable collection of items to be yielded.
    :param delay: delay (seconds) to wait until yielding the next item in the collection.
    :return:
    """
    for item in iterable:
        time.sleep(delay)
        yield item.proto


class ImdClient:
    """
    A simple IMD client, primarily for testing the IMD server.
    """

    def __init__(self, *, address: str, port: int):
        self.channel = grpc.insecure_channel("{0}:{1}".format(address, port))
        self.stub = InteractiveMolecularDynamicsStub(self.channel)
        self.threads = futures.ThreadPoolExecutor(max_workers=10)

    def publish_interactions_async(self, interactions) -> Future:
        """
        Publishes the collection of interactions on a thread, with an optional delay between
        each publication.
        :param interactions: A generator of interactions.
        """
        return self.threads.submit(self.publish_interactions, interactions)

    def publish_interactions(self, interactions) -> InteractionEndReply:
        """
        Publishes the collection of interactions on a thread, with an optional delay between
        each publication.
        :param interactions: A generator of interactions.
        :return: A reply indicating successful publishing of interaction.
        """
        return self.stub.PublishInteraction(interactions)

    def close(self):
        self.channel.close()
        self.threads.shutdown(wait=False)
