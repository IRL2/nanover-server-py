import asyncio
import json
import time
from functools import partial
from typing import Callable, Any, Optional

import aiohttp
import msgpack
import websockets
from nanover.app import NanoverImdApplication
from nanover.omni.converter import pack_array
from nanover.trajectory.frame_data import PARTICLE_COUNT, CHAIN_COUNT, RESIDUE_COUNT, SIMULATION_COUNTER, \
    PARTICLE_POSITIONS, PARTICLE_ELEMENTS, PARTICLE_RESIDUES, BOND_PAIRS, RESIDUE_CHAINS, BOX_VECTORS, FrameData
from nanover.utilities.cli import CancellationToken
from websockets import WebSocketServerProtocol


class NanoverServerWS:
    def __init__(
        self,
        *,
        app: NanoverImdApplication,
    ):
        self.app = app

    async def serve(
        self,
        *,
        name=None,
        ssl=None,
        cancellation: Optional[CancellationToken]=None,
    ):
        name = name or "nanover websocket server"

        async def await_cancellation():
            while not cancellation.is_cancelled:
                await asyncio.sleep(0.01)

        print("POLLING DISCOVERY")
        async with aiohttp.ClientSession() as session:
            response = await session.get("https://irl-discovery.onrender.com/list")
            print("list:", await response.json())

        print("RUNNING SERVER")
        async with websockets.serve(self.send_frames, "0.0.0.0", 0, ssl=ssl) as server:
            ip = get_local_ip()

            data = {
                "name": name,
                "web": f"https://{ip}:5500",
                "https": f"https://{ip}:5500",
            }

            if ssl is not None:
                port = list(server.sockets)[0].getsockname()[1]
                data["wss"] = f"wss://{ip}:{port}"
            else:
                port = list(server.sockets)[0].getsockname()[1]
                data["ws"] = f"ws://{ip}:{port}"

            async with websockets.connect("wss://irl-discovery.onrender.com/", open_timeout=5) as discovery:
                await asyncio.gather(
                    run_discovery(discovery, data),
                    await_cancellation(),
                )

    async def send_frames(self, websocket: WebSocketServerProtocol):
        print("INCOMING CONNECTION")

        async def send_frames():
            frame_publisher = self.app._frame_publisher

            prev = time.perf_counter()
            async for response in frame_publisher.subscribe_latest_frames():
                frame = FrameData(response.frame)
                data = {"frame": convert_frame(frame)}
                bytes = msgpack.packb(data)
                await websocket.send(bytes)
                next = time.perf_counter()
                print(next-prev, websocket)
                prev = next

        await send_frames()


async def run_discovery(websocket, data):
    init = json.loads(await websocket.recv())
    print(init, data)
    await websocket.send(json.dumps(data))


def get_local_ip():
    import socket;

    def attempt():
        yield from [
            ip
            for ip in socket.gethostbyname_ex(socket.gethostname())[2]
            if ip.startswith("192.")
        ][:1]

        s = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
        s.connect(("8.8.8.8", 53))
        name = s.getsockname()[0]
        s.close()
        yield name

    ip = next(attempt())

    return ip


pack_float32 = partial(pack_array, "f")
pack_uint32 = partial(pack_array, "I")
pack_uint8 = partial(pack_array, "B")

converters: dict[str, Callable[[], Any]] = {
    PARTICLE_COUNT: int,
    CHAIN_COUNT: int,
    RESIDUE_COUNT: int,
    SIMULATION_COUNTER: int,

    PARTICLE_POSITIONS: pack_float32,
    PARTICLE_ELEMENTS: pack_uint8,
    PARTICLE_RESIDUES: pack_uint32,

    BOND_PAIRS: pack_uint32,

    RESIDUE_CHAINS: pack_uint32,
    BOX_VECTORS: pack_float32,
}


def convert_frame(frame: FrameData):
    data = {}

    for key in frame.value_keys:
        converter = converters.get(key, lambda value: value)
        data[key] = converter(frame.values[key])

    for key in frame.array_keys:
        converter = converters.get(key, list)
        data[key] = converter(frame.arrays[key])

    return data
