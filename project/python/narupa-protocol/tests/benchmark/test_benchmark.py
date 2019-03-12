from narupy.benchmark.benchmark_client import BenchmarkClient
from narupy.benchmark.benchmark_server import BenchmarkServer, ServerCredentials
from datetime import datetime
import argparse
import numpy as np
import time
import os.path

path_to_creds = 'certification'
key_file = os.path.join(path_to_creds, '127.0.0.1.key')
cert_file = os.path.join(path_to_creds, '127.0.0.1.crt')
credentials = ServerCredentials(key_file, cert_file)

def test_insecure_server():
    host = "127.0.0.1:8001"
    server = BenchmarkServer(host)
    server.start()

def test_secure_server():
    host = "127.0.0.1:8002"
    server = BenchmarkServer(host, secure=True, credentials=credentials)
    server.start()

def test_insecure_connection():
    host = "127.0.0.1:8003"
    server = BenchmarkServer(host)
    server.start()
    client = BenchmarkClient(host)
    result = client.get_simple_message()
    assert result.payload == 3

def test_secure_connection():
    host = "127.0.0.1:8004"
    server = BenchmarkServer(host, secure=True, credentials=credentials)
    server.start()
    client = BenchmarkClient(host, secure=True, credentials=cert_file)
    result = client.get_simple_message()
    assert result.payload == 3

def setup_client_server(host):
    server = BenchmarkServer(host)
    server.start()
    client = BenchmarkClient(host)
    return client, server

def test_stream_frames():
    host = "127.0.0.1:8005"
    client, server =setup_client_server(host)
    n_atoms = 2
    n_frames = 100
    frames_rendezvous = client.get_frames(n_atoms,n_frames)
    count = 0
    for i, frame_response in enumerate(frames_rendezvous):
        count += 1
        pos = frame_response.frame.arrays['atom.position']
        floats = pos.float_values
        position_array = floats.values
        assert len(position_array) == n_atoms * 3
        assert all(value == i for value in position_array)
    assert count == n_frames

def test_stream_frames_raw():
    host = "127.0.0.1:8006"
    client, server =setup_client_server(host)
    n_atoms = 2
    n_frames = 100
    frames_rendezvous = client.get_frames_raw(n_atoms, n_frames)
    count = 0
    for i, frame_response in enumerate(frames_rendezvous):
        count += 1
        position_array = frame_response.frame
        assert len(position_array) == n_atoms * 3
        assert all(value == i for value in position_array)
    assert count == n_frames

def test_stream_frames_throttled():
    host = "127.0.0.1:8007"
    client, server =setup_client_server(host)
    n_atoms = 2
    n_frames = 100
    frames_rendezvous = client.get_frames_throttled(n_atoms, n_frames)
    count = 0
    for i, frame_response in enumerate(frames_rendezvous):
        count += 1
        pos = frame_response.frame.arrays['atom.position']
        floats = pos.float_values
        position_array = floats.values
        assert len(position_array) == n_atoms * 3
        assert all(value == i for value in position_array)
    assert count == n_frames