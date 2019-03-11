import grpc

import narupa.protocol.benchmark.benchmark_pb2_grpc as benchmark_pb2_grpc
import narupa.protocol.benchmark.benchmark_pb2 as benchmark_pb2
import narupa.protocol.instance.get_frame_pb2 as get_frame_pb2

class BenchmarkClient():
    channel: grpc.Channel
    stub: benchmark_pb2_grpc.StreamProviderStub

    def __init__(self, host:str, secure=False, credentials=None):
        if secure is False:
            self.channel = grpc.insecure_channel(host)

        self.stub = benchmark_pb2_grpc.StreamProviderStub(self.channel)


    def close(self):
        self.channel.close()

    def get_frame(self, n_atoms):
        req = benchmark_pb2.GetFrameRequest(n_atoms=n_atoms, n_frames=1)
        self.stub.GetFrames(req)

    def get_frames(self, n_atoms, n_frames = 100):
        req = benchmark_pb2.GetFrameRequest(n_atoms=n_atoms, n_frames=n_frames)
        return self.stub.GetFrames(req)

    def get_frames_gen(self, n_atoms, n_frames = 100):
        req = benchmark_pb2.GetFrameRequest(n_atoms=n_atoms, n_frames=n_frames)
        yield from self.stub.GetFrames(req)

    def get_simple_message(self):
        req = benchmark_pb2.SimpleMessage()
        req.payload = 3
        return self.stub.GetSimpleMessage(req)

    def get_frames_raw(self, n_atoms, n_frames = 100):
        req = benchmark_pb2.GetFrameRequest(n_atoms=n_atoms, n_frames=n_frames)
        return self.stub.GetFramesRaw(req)

    def get_frames_throttled(self,n_atoms, n_frames = 100):
        req = benchmark_pb2.GetFrameRequest(n_atoms=n_atoms, n_frames=n_frames)
        yield from self.stub.GetFramesThrottled(req)

