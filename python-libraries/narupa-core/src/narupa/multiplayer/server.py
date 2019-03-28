"""
Copyright (c) Mike O'Connor, University Of Bristol. All rights reserved.
Licensed under the GPL. See License.txt in the project root for license information.
"""

from concurrent import futures
from typing import List, Iterator
from narupa.protocol.multiplayer.multiplayer_pb2_grpc import MultiplayerServicer
from narupa.protocol.multiplayer.multiplayer_pb2 import PublishAvatarReply, Avatar, SceneProperties, ScenePropertyReply, ScenePropertyRequest
import narupa.protocol.multiplayer.multiplayer_pb2_grpc as mult_proto_grpc
import narupa.protocol.multiplayer.multiplayer_pb2 as mult_proto
from asyncio import Queue
import narupa.async.grpc_asyncio as grpc_asyncio
import grpc
import aiogrpc
from .multiplayer.scene import OwnerLock

class SceneProperties(object):
    """
    Represents properties of a multiplayer scene
    """
    _scene_property_lock : OwnerLock
    scene_properties : mult_proto.SceneProperties


    def __init__(self):
        super.__init__()
        self._scene_property_lock = OwnerLock()

    def is_locked(self):
        return self._scene_property_lock.is_locked()

    def try_lock(self, guid):
        return self._scene_property_lock.try_lock(guid)

    def try_unlock(self, guid):
        return self._scene_property_lock.try_unlock(guid)

    def is_lock_owner(self, guid):
        return self._scene_property_lock.is_locker(guid)


    def set_properties(self, properties:mult_proto.SceneProperties, lock_guid: int) -> bool:
        allow_set = True
        if self.is_locked():
            if not self.is_lock_owner(lock_guid):
                allow_set = False
        if allow_set:
            self.scene_properties = properties
            return True
        return False

class MultiplayerService(MultiplayerServicer):
    """
    Represents a multiplayer service.
    """
    avatar_queues: List[Queue]
    scene_property_queues: List[Queue]
    max_queue_size: int
    scene_properties: SceneProperties

    def __init__(self, max_avatar_queue_size:int =10):
        """
        Initialises the multiplayer service
        :param max_avatar_queue_size: Maximum number of avatar requests to place in queue.
        """
        super().__init__()

        self.avatar_queues = []
        self.max_queue_size = max_avatar_queue_size
        self.scene_properties = SceneProperties()


    async def SubscribeStream(self, request, context, queues: List[Queue]):
        queue = Queue()
        queues.append(queue)

        while context.is_active():
            # wait for the next avatar and return
            yield await queue.get()

    async def PublishStream(self, request_iterator: Iterator, context, queues):
        # loop over all requests
        # TODO async version of looping over iterator. this may block. https://blogs.gentoo.org/zmedico/2016/09/17/adapting-regular-iterators-to-asynchronous-iterators-in-python/
        for request in request_iterator:
            # loop over all subscribed queues
            for queue in queues:
                await queue.put(request)

    async def PublishMessage(self, message, context, queues):
        for queue in queues:
            await queue.put(message)

    async def SubscribeToAvatars(self, request, context) -> Avatar:
        """
        Subscribe to the avatar stream.
        :param request: Request to join avatar stream.
        :param context: 
        :return: Avatar
        """
        self.SubscribeStream(request, context, self.avatar_queues)

    async def PublishAvatar(self, request_iterator : Iterator, context) -> PublishAvatarReply:
        """
        Publish an avatar to the avatar stream
        :param request_iterator: Stream of avatars.
        :param context: 
        :return: 
        """

        self.PublishStream(request_iterator, context, self.avatar_queues)
        return PublishAvatarReply()

    async def SubscribeToSceneProperties(self, request, context) -> SceneProperties:
        """
        Subscribe to changes in the scene properties
        :param request: Request to join 
        :param context: 
        :return: 
        """
        self.SubscribeStream(self, request, context, self.scene_property_queues)

    async def SetLockScene(self, request: mult_proto.LockSceneProperty, context) -> mult_proto.LockSceneProperty:
        """
        Locks the scene for setting properties.
        A scene can only be locked by one user a time, providing their Guid for ownership. 
        :param request: Request to lock scene
        :param context: 
        :return: Lock scene property indicating whether scene was locked.
        """
        guid = request.guid
        is_locked = self.scene_properties.try_lock(guid)
        reply = mult_proto.LockSceneProperty()
        reply.locked = is_locked
        reply.guid = guid
        return reply

    async def SetSceneProperty(self, request : mult_proto.ScenePropertyRequest, context):
        """
        Sets the scene properties. 
        The scene can only be modified if it has been locked by the requester. 
        :param request: Request to set scene properties.
        :param context: 
        :return: Reply indicating whether scene property set was successful.
        """

        guid = request.property_lock_guid
        properties = request.properties
        success = self.scene_properties.set_properties(properties, guid)
        reply = ScenePropertyReply()
        reply.property_set = success
        reply.property_lock_guid = guid

        return reply

    async def UnlockScene(self, request, context):
        """
        Request to unlock the scene. 
        A scene can only be unlocked by whoever locked it, until the lock period has elapsed. 
        :param request: Request to unlock the scene.
        :param context: 
        :return: Reply indicating whether the scene property was unset. 
        """
        guid = request.guid
        is_unlocked = self.scene_properties.try_unlock(guid)
        reply = mult_proto.LockSceneProperty()
        reply.locked = not is_unlocked
        reply.guid = guid
        return reply

    async def JoinMultiplayer(self, request, context):
        # TODO implement
        pass

class MultiplayerServer(object):
    """
    Server providing multiplayer synchronisation.
    """

    def __init__(self, port='[::]:50051'):
        self.server = grpc.server(grpc_asyncio.AsyncioExecutor())
        self.multiplayer_services = MultiplayerService()
        mult_proto_grpc.add_MultiplayerServicer_to_server(
            self.multiplayer_services, self.server)
        self.server.add_insecure_port(port)
        self.server.start()

class MultiplayerClient(object):
    """
    Represents a client to the multiplayer server. Uses aiogrpc for asynchronous client behaviour.
    """
    stub:mult_proto_grpc.MultiplayerStub
    current_avatars:dict
    scene_properties:mult_proto.SceneProperties

    def __init__(self, channel:aiogrpc.Channel):
        self.stub = mult_proto_grpc.MultiplayerStub(channel)

    def join_multiplayer(self):
        # TODO implement
        pass

    async def join_avatar_stream(self):
        # TODO implement
        pass

    async def update_avatar(self, avatar):
        # TODO implement
        pass

    async def join_scene_properties_stream(self):
        # TODO implement
        pass

    async def set_scene_properties(self, properties):
        # TODO implement
        pass





