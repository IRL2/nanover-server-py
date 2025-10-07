import time
from contextlib import contextmanager
from threading import Lock

from nanover.utilities.queues import FrameMergingQueue
from nanover.utilities.timing import yield_interval
from nanover.trajectory import FrameData

SENTINEL = None


class FramePublisher:
    """
    An implementation of a trajectory service. Call send_frame
    to send data to clients when called by other python code.
    """

    def __init__(self):
        self.last_frame = FrameData()
        self.simulation_counter = 0

        self._frame_queues: set[FrameMergingQueue] = set()
        self._frame_queues_lock = Lock()
        self._last_frame_lock = Lock()
        self._request_id_lock = Lock()

    def subscribe_latest_frames(self, *, frame_interval=1 / 30, cancellation):
        """
        Yield the most recent frame, if changed, at a regular interval. Terminates when cancellation token is
        cancelled.
        """

        if cancellation.is_cancelled:
            return

        with self._last_frame_lock:
            initial_frame = self.last_frame

        with self._get_queue() as queue:
            if initial_frame is not None:
                yield initial_frame

            cancellation.subscribe_cancellation(lambda: queue.close())
            cancellation.subscribe_cancellation(lambda: print("CANCELLED"))

            for dt in yield_interval(frame_interval):
                item = queue.get(block=True)
                if cancellation.is_cancelled or item is SENTINEL:
                    break
                yield item

    @contextmanager
    def _get_queue(self):
        queue = FrameMergingQueue()

        try:
            with self._frame_queues_lock:
                self._frame_queues.add(queue)

            yield queue
        finally:
            queue.close()
            with self._frame_queues_lock:
                self._frame_queues.discard(queue)

    def _yield_last_frame_if_any(self):
        """
        Yields the last frame as a :class:`GetFrameResponse` object if there is
        one.

        This method places a lock on :attr:`last_frame` to prevent other threads modifying it as we
        read them.
        """
        with self._last_frame_lock:
            if self.last_frame is not None:
                yield self.last_frame

    def send_frame(self, frame_index: int, frame: FrameData):
        assert isinstance(frame, FrameData), "Frame must be of type FrameData"

        frame.server_timestamp = time.monotonic()
        frame.frame_index = frame_index

        if frame_index == 0:
            frame.simulation_counter = self.simulation_counter
            self.simulation_counter += 1

        with self._last_frame_lock:
            self.last_frame.update(frame)

        with self._frame_queues_lock:
            for queue in self._frame_queues:
                queue.put(frame)

    def close(self):
        with self._frame_queues_lock:
            for queue in self._frame_queues:
                queue.close()
