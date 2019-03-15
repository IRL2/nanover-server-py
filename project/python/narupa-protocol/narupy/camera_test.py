from narupa.protocol.renderer.camera_service_pb2 import CreateCameraRequest, UpdateCameraRequest
from narupa.protocol.renderer.camera_service_pb2_grpc import CameraServiceStub
from google.protobuf.json_format import MessageToJson, ParseDict
from narupa.protocol.instance.instance_service_pb2 import LoadTrajectoryRequest
from narupa.protocol.instance.instance_service_pb2_grpc import InstanceServiceStub

import grpc
import time
import random


def run():

    with grpc.insecure_channel('localhost:50051') as channel:
        camera = CameraServiceStub(channel)
        instance = InstanceServiceStub(channel)

        msg = LoadTrajectoryRequest()
        ParseDict({'instance_id': 'traj', 'path': 'https://files.rcsb.org/download/6AGY.pdb'}, msg)
        instance.LoadTrajectory(msg)


    def a():
        msg = CreateCameraRequest()
        ParseDict({'id': 'camera', 'properties': { 'background-color': '#FF0000'}}, msg)
        camera.CreateCamera(msg)

        while True:
            msg = UpdateCameraRequest()
            r = lambda: random.randint(0, 255)
            color = '#%02X%02X%02X' % (r(), r(), r())
            ParseDict({'id': 'camera', 'properties': {'background-color': color}}, msg)
            camera.UpdateCamera(msg)

            time.sleep(0.1)




run()