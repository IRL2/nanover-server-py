import grpc
from concurrent import futures
import narupa.protocol.instance.molecule_provider_pb2_grpc as foo


class GrpcServer:

    def __init__(self, *, address : str, port : int):
        self.server = grpc.server(futures.ThreadPoolExecutor(max_workers=10))

        self.setup_services()

        self.server.add_insecure_port(address="{0}:{1}".format(address, port))
        self.server.start()

    def setup_services(self):
        pass
