import argparse
from datetime import datetime
import os

import grpc

import narupa.protocol.benchmark.benchmark_pb2_grpc as benchmark_pb2_grpc
import narupa.protocol.benchmark.benchmark_pb2 as benchmark_pb2
import narupa.protocol.instance.get_frame_pb2 as get_frame_pb2
import numpy as np

class BenchmarkClient():
    channel: grpc.Channel
    stub: benchmark_pb2_grpc.StreamProviderStub

    def __init__(self, host:str, secure=False, credentials=None):
        if secure is False:
            self.channel = grpc.insecure_channel(host)
        else:
            with open(credentials, 'rb') as file:
                creds = grpc.ssl_channel_credentials(file.read())
            self.channel = grpc.secure_channel(host, creds)
        self.stub = benchmark_pb2_grpc.StreamProviderStub(self.channel)


    def close(self):
        self.channel.close()

    def get_frame(self, n_atoms):
        req = benchmark_pb2.GetFrameRequest(n_atoms=n_atoms, n_frames=1)
        self.stub.GetFrames(req)

    def get_frames(self, n_atoms, n_frames = 100):
        req = benchmark_pb2.GetFrameRequest(n_atoms=n_atoms, n_frames=n_frames)
        return self.stub.GetFrames(req)

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

def time_frames(method, n_atoms, n_frames):
    times = []
    current_time = datetime.now()
    frames_received = 0
    for frame in method(n_atoms, n_frames):
        previous_time, current_time = current_time, datetime.now()
        elapsed_ms = (current_time - previous_time).total_seconds() * 1000
        times.append(elapsed_ms)
        frames_received += 1
    print('method', method, 'received: ', frames_received, ' avg: ', np.mean(times), ' sd: ', np.std(times))
    return frames_received, times

def run(args):
    client = BenchmarkClient(args.host, secure=args.secure, credentials=args.server_certificate_file)

    time_frames(client.get_frames, args.n_atoms, args.n_frames)
    time_frames(client.get_frames_throttled, args.n_atoms, args.n_frames)
    time_frames(client.get_frames_raw, args.n_atoms, args.n_frames)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Benchmark a gRPC client receiving atomic data from a server')
    parser.add_argument('n_atoms', action='store_const', help='Number of atoms to simulate data', const=2048)
    parser.add_argument('n_frames', action='store_const', help='Number of frames to run.', const=500)
    parser.add_argument('host', action='store_const', help='Host address', const='127.0.0.1:8007')
    path_to_creds = '../../../../../certification'
    parser.add_argument('secure', action='store_const', help='Whether to run securely', const=True)
    parser.add_argument('server_certificate_file', action='store_const', help='Server certificate file', const=os.path.join(path_to_creds, '127.0.0.1.crt'))

    args = parser.parse_args()
    run(args)
