import argparse
from datetime import datetime, timedelta
import time
from concurrent import futures

import grpc
import numpy as np

import narupa.protocol.benchmark.benchmark_pb2_grpc as benchmark_pb2_grpc
from narupa.protocol.benchmark.benchmark_pb2 import RawFrame
import narupa.protocol.instance.get_frame_pb2 as get_frame_pb2
import narupa.protocol.trajectory.frame_pb2 as frame_pb2


class BenchmarkService(benchmark_pb2_grpc.StreamProviderServicer):
    n_frames: int
    frame_request: get_frame_pb2.GetFrameRequest


    def __init__(self, fps=30):
        self.array = np.zeros(2048 * 3, dtype=np.float32)
        self.array.fill(1.0)
        self.frame = frame_pb2.FrameData()
        self.frame_request = get_frame_pb2.GetFrameResponse(frame_index=0, frame=self.frame)
        # this is really gross and needs wrapping??
        self.valuearray = self.frame_request.frame.arrays['atom.position'].float_values.values
        self.valuearray.extend(self.array)
        self.raw_frame = RawFrame()
        self.raw_frame.frame.extend(self.array)
        self.fps = 30



    def GetSimpleMessage(self, request, context):
        return request

    def FillFrame(self, array):
        self.valuearray[:] = array
        return self.frame_request

    def ServeFrame(self, array, value):
        self.array.fill(value)
        return self.FillFrame(self.array)



    def GetFrames(self, request, context):
        print(request)
        self.array.resize(request.n_atoms * 3)
        for i in range(request.n_frames):
            self.array.fill(i)
            yield self.FillFrame(self.array)

    def GetFramesRaw(self, request, context):
        self.array.resize(request.n_atoms * 3)
        for i in range(request.n_frames):
            self.array.fill(i)
            self.raw_frame.frame[:] = self.array
            yield self.raw_frame

    def GetFramesThrottled(self, request, context):
        self.array.resize(request.n_atoms * 3)
        target_delta = timedelta(milliseconds= 1.0 / self.fps * 1000)
        current_time = target_time = datetime.now()
        target_time += target_delta

        for i in range(request.n_frames):
            # compute frame
            frame = self.ServeFrame(self.array, i)
            # determine how much time has elapsed since last frame
            previous_time, current_time = current_time, datetime.now()

            sleep_time = target_time - current_time
            sleep_seconds = sleep_time.total_seconds()
            if sleep_seconds > 0:
                print('sleeping for ', sleep_seconds)
                time.sleep(sleep_time.seconds)

            current_time = datetime.now()
            target_time = current_time + target_delta
            yield frame



class BenchmarkServer():
    server: grpc.Server

    def __init__(self, host: str, secure: bool=False, credentials=None):
        self.server = grpc.server(futures.ThreadPoolExecutor(max_workers=10))
        benchmark_pb2_grpc.add_StreamProviderServicer_to_server(BenchmarkService(), self.server)

        if not secure:
            self.server.add_insecure_port(host)
        else:
            pass

    def start(self):
        self.server.start()

    def stop(self):
        self.server.stop(0)

def run(args):
    server = BenchmarkServer(args.host, secure=False)
    server.start()
    try:
        while True:
            time.sleep(10000)
    except KeyboardInterrupt:
        server.stop()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Benchmark a gRPC client and server running over localhost.')
    parser.add_argument('host', action='store_const', help='Host address', const='[::]:8001')
    args = parser.parse_args()
    run(args)