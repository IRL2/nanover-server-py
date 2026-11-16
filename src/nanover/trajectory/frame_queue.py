from nanover.utilities.queues import LastItemQueue

from .frame_wrapper import FrameData


class FrameMergingQueue(LastItemQueue[FrameData]):
    """
    SingleItemQueue specifically for FrameData items. Put frames will be
    aggregated with any existing frame so that there is at most one frame in the
    queue at any time.
    """

    def close(self):
        self.put(None)

    def put(self, item: FrameData | None, **kwargs):
        with self._lock:
            if item is None:
                # The None sentinel value indicates that the queue user should terminate,
                # so it is safe to discard aggregated frames.
                self._item = None
            else:
                if self._item is None:
                    self._item = FrameData()
                self._item.update(item)
            self._has_item = True
            self.not_empty.notify()
