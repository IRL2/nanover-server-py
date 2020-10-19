# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.
import time
from queue import Queue
from threading import Lock
from typing import Union, Callable

from narupa.utilities.request_queues import (
    DictOfQueues, GetFrameResponseAggregatingQueue)
from narupa.utilities.timing import yield_interval
from narupa.protocol.trajectory import (
    TrajectoryServiceServicer, GetFrameResponse,
    add_TrajectoryServiceServicer_to_server)
from narupa.protocol.trajectory import FrameData as RawFrameData
from narupa.trajectory.frame_data import FrameData, SERVER_TIMESTAMP

SENTINEL = None

FRAME_SERVICE_NAME = "trajectory"


class FramePublisher(TrajectoryServiceServicer):
    """
    An implementation of a trajectory service. Call send_frame
    to send data to clients when called by other python code.
    """

    frame_queues: DictOfQueues
    last_frame: RawFrameData
    last_frame_index: int
    last_request_id: int
    _frame_queue_lock: Lock
    _last_frame_lock: Lock
    _request_id_lock: Lock

    def __init__(self):
        self.name: str = FRAME_SERVICE_NAME
        self.add_to_server_method: Callable = add_TrajectoryServiceServicer_to_server
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
        yield from self._subscribe_frame_base(request,
                                              context,
                                              queue_type=Queue)

    def SubscribeLatestFrames(self, request, context):
        """
        Subscribe to the last produced frames produced by the service.

        This method publishes the latest frame available at the time of
        yielding.
        """
        yield from self._subscribe_frame_base(
            request,
            context,
            queue_type=GetFrameResponseAggregatingQueue,
        )

    def _subscribe_frame_base(self, request, context, queue_type):
        listen_for_cancellation = context.add_callback
        request_id = self._get_new_request_id()
        yield from self._yield_last_frame_if_any()
        with self.frame_queues.one_queue(request_id, queue_class=queue_type) as queue:
            if not listen_for_cancellation(lambda: queue.put(SENTINEL)):
                return
            for dt in yield_interval(request.frame_interval):
                item = queue.get(block=True)
                if item is SENTINEL:
                    break
                yield item

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
        Yields the last frame as a :class:`GetFrameResponse` object if there is
        one.

        This method places a lock on :attr:`last_frame` and
        :attr:`last_frame_index` to prevent other threads to modify them as we
        read them.
        """
        with self._last_frame_lock:
            if self.last_frame is not None:
                yield GetFrameResponse(frame_index=self.last_frame_index, frame=self.last_frame)

    def send_frame(self, frame_index: int, frame: Union[FrameData, RawFrameData]):
        now = time.monotonic()
        if isinstance(frame, FrameData):
            frame.server_timestamp = now
            frame = frame.raw
        else:
            frame.values[SERVER_TIMESTAMP].number_value = now

        with self._last_frame_lock:
            if self.last_frame is None:
                self.last_frame = RawFrameData()
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
            queue.put(SENTINEL)
