from __future__ import print_function
import logging

import grpc

import narupa.protocol.multiplayer.async_test_pb2_grpc as async_grpc
import narupa.protocol.multiplayer.async_test_pb2 as async_test


def run():
    # NOTE(gRPC Python Team): .close() is possible on a channel and should be
    # used in circumstances in which the with statement does not fit the needs
    # of the code.
    with grpc.insecure_channel('localhost:9876') as channel:
        stub = async_grpc.TestStub(channel)
        response = stub.DoSomething(async_test.Request())

        for r in response:
            print('response:', r)


if __name__ == '__main__':
    logging.basicConfig()
    run()
