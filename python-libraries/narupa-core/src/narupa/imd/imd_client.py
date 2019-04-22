import time
from concurrent import futures

import grpc

from narupa.protocol.imd import InteractiveMolecularDynamicsStub


class ImdClient:

    def __init__(self, *, address: str, port: int):
        self.channel = grpc.insecure_channel("{0}:{1}".format(address, port))
        self.stub = InteractiveMolecularDynamicsStub(self.channel)
        self.threads = futures.ThreadPoolExecutor(max_workers=10)

    def _to_generator(self, list, delay=0):
        for item in list:
            time.sleep(delay)
            yield item.proto

    def publish_interactions_async(self, interactions, delay):
        self.threads.submit(self.publish_interactions, interactions, delay)

    def publish_interactions(self, interactions, delay):
        return self.stub.PublishInteraction(self._to_generator(interactions,delay))

    def close(self):
        self.channel.close()
        self.threads.shutdown(wait=False)
