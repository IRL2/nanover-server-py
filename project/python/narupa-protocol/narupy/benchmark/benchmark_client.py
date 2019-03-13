"""
Copyright (c) Mike O'Connor, University Of Bristol. All rights reserved.
Licensed under the GPL. See License.txt in the project root for license information.

Benchmark client for narupa style data.
"""

import argparse
import csv
from collections import namedtuple
from datetime import datetime
import os

import grpc

import narupa.protocol.benchmark.benchmark_pb2_grpc as benchmark_pb2_grpc
import narupa.protocol.benchmark.benchmark_pb2 as benchmark_pb2
import narupa.protocol.instance.get_frame_pb2 as get_frame_pb2
import numpy as np

BenchmarkRun = namedtuple('BenchmarkRun', 'client, n_atoms, n_frames, secure, output')

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
    print('method', method.__name__, 'received: ', frames_received, ' avg: ', np.mean(times), ' sd: ', np.std(times))
    return frames_received, times

def write_csv_header(csv_path):
    with open(csv_path, 'w') as f:
        writer = csv.writer(f)
        header = ['time', 'method', 'n_atoms', 'n_frames', 'mean', 'std', 'secure']
        writer.writerow(header)


def run_benchmark(args : BenchmarkRun):
    client = args.client
    methods = [client.get_frames, client.get_frames_throttled, client.get_frames_raw]
    with open(args.output, 'a') as f:
        writer = csv.writer(f)
        for method in methods:
            frames, times = time_frames(method, args.n_atoms, args.n_frames)
            row = [datetime.now(), method.__name__, args.n_atoms, args.n_frames, np.mean(times), np.std(times), args.secure]
            writer.writerow(row)

def run(args):
    client = BenchmarkClient(args.host, secure=args.secure, credentials=args.server_certificate_file)

    if not os.path.exists(args.output):
        write_csv_header(args.output)
    range = [int(x) for x in np.logspace(args.n_atoms_min, args.n_atoms_max, args.n_samples, base=2)]
    print(range)
    for n_atoms in range:
        benchmark_args = BenchmarkRun(client=client, n_atoms=n_atoms, n_frames=args.n_frames, output=args.output, secure=args.secure)
        print('running: ', benchmark_args)
        run_benchmark(benchmark_args)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Benchmark a gRPC client receiving atomic data from a server')
    parser.add_argument('n_atoms_min', action='store_const', help='Number of atoms to simulate data (base 2)', const=10)
    parser.add_argument('n_atoms_max', action='store_const', help='Number of atoms to simulate data (base 2)', const=15)
    parser.add_argument('n_samples', action='store_const', help='Number of samples in atom stride', const=8)
    parser.add_argument('n_frames', action='store_const', help='Number of frames to run.', const=5000)
    parser.add_argument('host', action='store_const', help='Host address', const='0.0.0.0:8000')
    path_to_creds = '../../../../../certification'
    parser.add_argument('secure', action='store_const', help='Whether to run securely', const=False)
    parser.add_argument('server_certificate_file', action='store_const', help='Server certificate file', const=os.path.join(path_to_creds, '127.0.0.1.crt'))
    parser.add_argument('output', action='store_const', const="results.csv")
    args = parser.parse_args()
    run(args)
