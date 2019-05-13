from queue import Queue, Empty
from typing import List

from narupa.protocol.trajectory import TrajectoryServiceServicer, GetFrameResponse, FrameData


class FramePublisher(TrajectoryServiceServicer):
    """
    An implementation of a trajectory service. Call send_frame
    to send data to clients when called by other python code.
    """

    frame_queues: List[Queue]

    last_frame: FrameData

    def __init__(self):
        self.frame_queues = []
        self.last_frame = None
        self.last_frame_index = 0

    def SubscribeFrames(self, request, context):

        if self.last_frame is not None:
            yield GetFrameResponse(frame_index=self.last_frame_index, frame=self.last_frame)

        queue = Queue()
        self.frame_queues.append(queue)

        while True:
            try:
                item = queue.get(block=True, timeout=0.5)
            except Empty:
                pass
            else:
                yield item
            finally:
                queue.task_done()

    def send_frame(self, frame_index: int, frame: FrameData):
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

        for queue in self.frame_queues:
            queue.put(GetFrameResponse(frame_index=frame_index, frame=frame))
