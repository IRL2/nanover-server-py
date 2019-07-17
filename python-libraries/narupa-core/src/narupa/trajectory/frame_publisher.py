from queue import Queue, Empty
from threading import Lock

from narupa.core.request_queues import DictOfQueues, SingleItemQueue
from narupa.protocol.trajectory import TrajectoryServiceServicer, GetFrameResponse, FrameData

import time


class FramePublisher(TrajectoryServiceServicer):
    """
    An implementation of a trajectory service. Call send_frame
    to send data to clients when called by other python code.
    """

    frame_queues: DictOfQueues
    last_frame: FrameData
    last_frame_index: int
    last_request_id: int
    _frame_queue_lock: Lock
    _last_frame_lock: Lock
    _request_id_lock: Lock

    def __init__(self):
        self.frame_queues = DictOfQueues()
        self.last_frame = None
        self.last_frame_index = 0
        self.last_request_id = 0
        self._last_frame_lock = Lock()
        self._request_id_lock = Lock()

    def SubscribeFrames(self, request, context):
        """
        Subscribe to all the frames produced by the service.

        This method publishes all the frames produced by the trajectory service,
        starting when the client subscribes.
        """
        yield from self._subscribe_frame_base(request, context, queue_type=Queue)

    def SubscribeLatestFrames(self, request, context):
        """
        Subscribe to the last produced frames produced by the service.

        This method publishes the latest frame available at the time of yielding.
        """
        yield from self._subscribe_frame_base(request,
                                              context,
                                              queue_type=SingleItemQueue,
                                              frame_delay=1/30)

    def _subscribe_frame_base(self, request, context, queue_type, frame_delay=None):
        request_id = self._get_new_request_id()
        yield from self._yield_last_frame_if_any()

        with self.frame_queues.one_queue(request_id, queue_class=queue_type) as queue:
            while context.is_active():
                item = queue.get(block=True)
                if item is None:
                    break
                yield item
                if frame_delay:
                    time.sleep(frame_delay)

    def _get_new_request_id(self) -> int:
        """
        Provide a new client id in a thread safe way.
        """
        with self._request_id_lock:
            self.last_request_id += 1
            client_id = self.last_request_id
        return client_id

    def _yield_last_frame_if_any(self):
        """
        Yields the last frame as a :class:`GetFrameResponse` object if there is one.

        This method places a lock on :attr:`last_frame` and
        :attr:`last_frame_index` to prevent other threads to modify them as we
        read them.
        """
        with self._last_frame_lock:
            if self.last_frame is not None:
                yield GetFrameResponse(frame_index=self.last_frame_index, frame=self.last_frame)

    def send_frame(self, frame_index: int, frame: FrameData):
        with self._last_frame_lock:
            if self.last_frame is None:
                self.last_frame = FrameData()
            self.last_frame_index = frame_index

            for key in frame.arrays.keys():
                if key in self.last_frame.arrays:
                    del self.last_frame.arrays[key]
            for key in frame.values.keys():
                if key in self.last_frame.values:
                    del self.last_frame.values[key]

            self.last_frame.MergeFrom(frame)

        for queue in self.frame_queues.iter_queues():
            queue.put(GetFrameResponse(frame_index=frame_index, frame=frame))

    def close(self):
        for queue in self.frame_queues.iter_queues():
            queue.put(None)
