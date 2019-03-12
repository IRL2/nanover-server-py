"""
Benchmark program for Narupa style data.
"""
from narupy.benchmark.benchmark_client import BenchmarkClient
from narupy.benchmark.benchmark_server import BenchmarkServer, ServerCredentials
from datetime import datetime
import argparse
import numpy as np
import os.path


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

def get_all_frames_unary(client: BenchmarkClient, n_atoms, n_frames):
    for i in range(n_frames):
        client.get_frame(n_atoms)


def get_all_frames_stream(client: BenchmarkClient, n_atoms, n_frames):
    return sum([1 for frame in client.get_frames(n_atoms, n_frames)])


def get_all_frames_stream_raw(client: BenchmarkClient, n_atoms, n_frames):
    return sum([1 for frame in client.get_frames_raw(n_atoms, n_frames)])


def get_simple_message(client: BenchmarkClient, repeats):
    for i in range(repeats):
        client.get_simple_message()

def get_all_frames_throttled(client: BenchmarkClient, n_atoms, n_frames):
    times = []
    current_time = datetime.now()
    frames_received = 0
    for frame in client.get_frames_throttled(n_atoms, n_frames):
        previous_time, current_time = current_time, datetime.now()
        elapsed_ms = (current_time - previous_time).total_seconds() * 1000
        times.append(elapsed_ms)
        frames_received += 1
    return frames_received, times


def run(args):
    creds = ServerCredentials(args.server_private_key, args.server_certificate_file)
    server = BenchmarkServer(args.host, secure=args.secure, credentials=creds)
    server.start()
    client = BenchmarkClient(args.host, secure=args.secure, credentials=args.server_certificate_file)

    time_frames(client.get_frames, args.n_atoms, args.n_frames)
    time_frames(client.get_frames_throttled, args.n_atoms, args.n_frames)
    time_frames(client.get_frames_raw, args.n_atoms, args.n_frames)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Benchmark a gRPC client and server running over localhost.')
    parser.add_argument('n_atoms', action='store_const', help='Number of atoms to simulate data', const=32000)
    parser.add_argument('n_frames', action='store_const', help='Number of frames to run.', const=200)
    parser.add_argument('host', action='store_const', help='Host address', const='127.0.0.1:8006')
    path_to_creds = '../../../../../certification'
    parser.add_argument('secure', action='store_const', help='Whether to run securely', const=True)
    parser.add_argument('server_private_key', action='store_const', help='Server private key file', const=os.path.join(path_to_creds, '127.0.0.1.key'))
    parser.add_argument('server_certificate_file', action='store_const', help='Server certificate file', const=os.path.join(path_to_creds, '127.0.0.1.crt'))

    args = parser.parse_args()
    run(args)
