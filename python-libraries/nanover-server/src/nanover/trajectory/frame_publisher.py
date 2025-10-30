import time
from contextlib import contextmanager
from threading import Lock

from nanover.utilities.cli import CancellationToken
from nanover.utilities.queues import FrameMergingQueue
from nanover.utilities.timing import yield_interval
from nanover.trajectory import FrameData


class FramePublisher:
    """
    Manages a single FrameData and provides a means to update it and make multiple independent subscriptions to a
    stream of those updates.
    """

    def __init__(self):
        self.simulation_counter = 0
        self.next_frame_index = 0

        self.current_frame: FrameData | None = None
        self._frame_queues: set[FrameMergingQueue] = set()
        self._frame_queues_lock = Lock()

    def close(self):
        """
        Terminate all frame subscriptions.
        """
        with self._frame_queues_lock:
            for queue in self._frame_queues:
                queue.close()

    def send_clear(self):
        """
        Queue the start of a new sequence of frames that don't reuse the topology etc of the previous frames.
        """
        self.next_frame_index = 1

        frame = FrameData()
        frame.frame_index = 0
        frame.simulation_counter = self.simulation_counter
        self.simulation_counter += 1

        self.current_frame = frame
        self._queue_frame_for_all_subscribers(frame)

    def send_frame(self, frame: FrameData):
        """
        Queue a frame for publishing.
        """
        if self.current_frame is None:
            self.send_clear()

        prev_frame = self.current_frame or FrameData()

        frame.server_timestamp = time.monotonic()
        frame.frame_index = self.next_frame_index
        self.next_frame_index += 1

        next_frame = FrameData()
        next_frame.update(prev_frame)
        next_frame.update(frame)

        self.current_frame = next_frame
        self._queue_frame_for_all_subscribers(frame)

    def subscribe_latest_frames(
        self, *, frame_interval=1 / 30, cancellation: CancellationToken
    ):
        """
        Yield the most recent frame, if changed, at a regular interval. Terminates when cancellation token is
        cancelled or this FramePublisher is closed.
        """
        if cancellation.is_cancelled:
            return

        with self._get_queue() as queue:
            cancellation.subscribe_cancellation(lambda: queue.close())

            for _ in yield_interval(frame_interval):
                item = queue.get(block=True)
                if cancellation.is_cancelled or item is None:
                    break
                yield item

    def _queue_frame_for_all_subscribers(self, frame):
        with self._frame_queues_lock:
            for queue in self._frame_queues:
                queue.put(frame)

    @contextmanager
    def _get_queue(self):
        queue = FrameMergingQueue()

        # insert the current aggregate frame if it exists
        if self.current_frame is not None:
            queue.put(self.current_frame)

        try:
            with self._frame_queues_lock:
                self._frame_queues.add(queue)
            yield queue
        finally:
            with self._frame_queues_lock:
                self._frame_queues.discard(queue)
            queue.close()
