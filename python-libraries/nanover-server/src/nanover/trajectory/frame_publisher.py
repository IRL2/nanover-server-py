import time
from threading import Lock

from nanover.utilities.request_queues import (
    DictOfQueues,
    GetFrameResponseAggregatingQueue,
)
from nanover.utilities.timing import yield_interval
from nanover.protocol.trajectory import FrameData as RawFrameData
from nanover.trajectory import FrameData, FrameData2
from nanover.trajectory.convert import convert_grpc_frame_to_dict_frame

SENTINEL = None


class FramePublisher:
    """
    An implementation of a trajectory service. Call send_frame
    to send data to clients when called by other python code.
    """

    last_frame_index: int
    last_request_id: int
    _frame_queue_lock: Lock
    _last_frame_lock: Lock
    _request_id_lock: Lock

    def __init__(self):
        self.frame_queues = DictOfQueues()
        self.last_frame = FrameData2()
        self.last_frame_index = 0
        self.last_request_id = 0
        self.simulation_counter = 0
        self._last_frame_lock = Lock()
        self._request_id_lock = Lock()

    def subscribe_latest_frames(
        self,
        *,
        frame_interval=1 / 30,
        cancellation,
        queue_class=GetFrameResponseAggregatingQueue
    ):
        """
        Yield the most recent frame, if changed, at a regular interval. Terminates when cancellation token is
        cancelled.
        """
        request_id = self._get_new_request_id()

        with self.frame_queues.one_queue(
            request_id,
            queue_class=queue_class,
        ) as queue:
            if cancellation.is_cancelled:
                return

            with self._last_frame_lock:
                initial_frame = self.last_frame

            if initial_frame is not None:
                yield initial_frame

            cancellation.subscribe_cancellation(lambda: queue.put(SENTINEL))

            for dt in yield_interval(frame_interval):
                item = queue.get(block=True)
                if cancellation.is_cancelled or item is SENTINEL:
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
                yield self.last_frame

    def send_frame(
        self, frame_index: int, frame: FrameData | RawFrameData | FrameData2
    ):
        actual_frame: FrameData2

        if isinstance(frame, RawFrameData):
            actual_frame = FrameData2(convert_grpc_frame_to_dict_frame(frame))
        elif isinstance(frame, FrameData):
            actual_frame = FrameData2(convert_grpc_frame_to_dict_frame(frame.raw))
        elif isinstance(frame, FrameData2):
            actual_frame = frame
        else:
            raise TypeError("Invalid frame type")

        actual_frame.server_timestamp = time.monotonic()
        actual_frame.frame_index = frame_index

        if frame_index == 0:
            actual_frame.simulation_counter = self.simulation_counter

        with self._last_frame_lock:
            self.last_frame_index = frame_index
            self.last_frame.update(actual_frame)

        for queue in self.frame_queues.iter_queues():
            queue.put(frame)

    def close(self):
        for queue in self.frame_queues.iter_queues():
            queue.put(SENTINEL)
