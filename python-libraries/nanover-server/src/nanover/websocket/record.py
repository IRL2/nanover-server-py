from concurrent.futures import ThreadPoolExecutor
from contextlib import suppress
from time import perf_counter_ns

import msgpack
from websockets import ConnectionClosed
from websockets.sync.client import connect
from websockets.sync.connection import Connection

from nanover.omni import OmniRunner
from nanover.recording.reading import MessageEvent
from nanover.recording.writing import NanoverRecordingWriter

from .client.app_client import get_websocket_address_from_app_server
from .client.base_client import MAX_MESSAGE_SIZE


class BackgroundRecordingContext:
    @classmethod
    def from_address_to_path(cls, *, address: str, path: str):
        """
        Connect to the given websocket address and record trajectory frames and state updates to a file at
        the given path

        :param address: Websocket address in the form protocol://host[:port]
        :param path: File path to record to
        """
        writer = NanoverRecordingWriter(path)
        connection = connect(address, max_size=MAX_MESSAGE_SIZE)
        return cls(connection, writer)

    def __init__(self, connection: Connection, writer: NanoverRecordingWriter):
        self._connection = connection
        self._writer = writer
        self._threads = ThreadPoolExecutor(max_workers=1)
        self._open = True

        self.future = self._threads.submit(self._record)
        self.future.add_done_callback(lambda _: self.close)

    def _record(self):
        with self._writer:
            for event in message_events_from_websocket(self._connection):
                self._writer.write_message_event(event)

    def close(self):
        """
        End recording by disconnection and close writer.
        """
        if not self._open:
            return

        self._open = False
        self._connection.close()
        self._writer.close()
        self._threads.shutdown()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()


def record_from_runner(runner: OmniRunner, out_path):
    """
    Connect to the given runner and record trajectory frames and state updates to a file

    :param runner: OmniRunner instance to connect to
    :param out_path: File to write recording to
    """
    return BackgroundRecordingContext.from_address_to_path(
        address=get_websocket_address_from_app_server(runner.app_server),
        path=out_path,
    )


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
