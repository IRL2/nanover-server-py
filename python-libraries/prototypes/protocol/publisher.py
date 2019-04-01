from __future__ import annotations
from typing import List
from queue import Queue, Empty

#
#
#
#

class Subscription(object):
    topics: List[str]
    publisher: Publisher

    message_queue: Queue[Message]

    def __init__(self, publisher: Publisher, id):
        self.publisher = publisher
        self.id = id
        self.message_queue = []


class Publisher(object):
    topics: List[str]
    subscriptions: List[Subscription]

    def __init__(self):
        pass

    def add_subscription(self, id):
        subscription: Subscription = Subscription(id)
        self.subscriptions.append(subscription)
        return subscription

    def publish(self, message):
        for subscription in self.subscriptions:
            subscription.message_queue.put(message)

    def start_stream(self, request, context):
        subscription = Subscription(context, request.topics)
        while True:
            item = subscription.message_queue.get(True)
            for packet in item.packets:
                yield packet

