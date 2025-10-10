from concurrent.futures import ThreadPoolExecutor
from contextlib import suppress
from time import perf_counter_ns

import msgpack
from websockets import ConnectionClosed
from websockets.sync.client import connect
from websockets.sync.connection import Connection

from nanover.omni import OmniRunner
from nanover.recording.reading import MessageEvent
from nanover.recording.writing import record_messages
from nanover.websocket.client.app_client import get_websocket_address_from_app_server
from nanover.websocket.client.base_client import MAX_MESSAGE_SIZE


def record_from_runner(runner: OmniRunner, out_path):
    """
    Connect to the given runner and record trajectory frames and state updates to a file

    :param runner: OmniRunner instance to connect to
    :param out_path: File to write recording to
    """
    return record_from_server(
        address=get_websocket_address_from_app_server(runner.app_server),
        out_path=out_path,
    )


def record_from_server(address, out_path):
    """
    Connect to the given address and record trajectory frames and state updates to a file

    :param address: Address (protocol://host:port) of server to connect to
    :param out_path: File to write recording to
    """

    websocket = connect(address, max_size=MAX_MESSAGE_SIZE)
    executor = ThreadPoolExecutor(max_workers=1)
    executor.submit(record_messages, out_path, message_events_from_websocket(websocket))

    return executor, websocket


def message_events_from_websocket(websocket: Connection):
    """
    Iterate the stream of incoming messages of a websocket connection and yield generic MessageEvents for use with the
    recording functions.
    """
    start_time = perf_counter_ns()

    def get_timestamp():
        return int((perf_counter_ns() - start_time) / 1000)

    with suppress(ConnectionClosed), websocket:
        for data in websocket:
            if isinstance(data, bytes):
                yield MessageEvent(
                    timestamp=get_timestamp(),
                    message=msgpack.unpackb(data),
                )
