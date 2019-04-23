# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.

import time
from concurrent import futures
from typing import Collection

import grpc

from narupa.imd.interaction import Interaction
from narupa.protocol.imd import InteractiveMolecularDynamicsStub


class ImdClient:
    """
    A simple IMD client, primarily for testing the IMD server.
    """

    def __init__(self, *, address: str, port: int):
        self.channel = grpc.insecure_channel("{0}:{1}".format(address, port))
        self.stub = InteractiveMolecularDynamicsStub(self.channel)
        self.threads = futures.ThreadPoolExecutor(max_workers=10)

    def _to_generator(self, list, delay: float=0):
        for item in list:
            time.sleep(delay)
            yield item.proto

    def publish_interactions_async(self, interactions: Collection[Interaction], delay: float=0):
        """
        Publishes the collection of interactions on a thread, with an optional delay between
        each publication.
        :param interactions: Collection of interactions.
        :param delay: Time to delay between publishing interactions, in seconds.
        :return:
        """
        self.threads.submit(self.publish_interactions, interactions, delay)

    def publish_interactions(self, interactions: Collection[Interaction], delay: float=0):
        """
        Publishes the collection of interactions on a thread, with an optional delay between
        each publication.
        :param interactions: Collection of interactions.
        :param delay: Time to delay between publishing interactions, in seconds.
        :return:
        """
        return self.stub.PublishInteraction(self._to_generator(interactions, delay))

    def close(self):
        self.channel.close()
        self.threads.shutdown(wait=False)
