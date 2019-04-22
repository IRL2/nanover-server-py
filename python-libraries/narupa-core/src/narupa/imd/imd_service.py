from queue import Queue
from typing import List

from narupa.protocol.imd import InteractiveMolecularDynamicsServicer


class ImdService(InteractiveMolecularDynamicsServicer):

    interaction_queues: List[Queue]

    def __init__(self):
        self.interaction_queues = []