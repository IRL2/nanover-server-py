from os import PathLike
from nanover.recording2.record2 import MessageZipReader
from nanover.websocket.convert import (
    unpack_dict_frame,
    convert_dict_frame_to_grpc_frame,
    convert_dict_state_to_dictionary_change,
)


def iter_recording_file(path: PathLike[str]):
    reader = MessageZipReader.from_path(path)
    for entry in reader:
        frame, update = None, None
        message = reader.get_message_from_entry(entry)
        if "frame" in message:
            frame = convert_dict_frame_to_grpc_frame(
                unpack_dict_frame(message["frame"])
            )
        if "state" in message:
            update = convert_dict_state_to_dictionary_change(message["state"])

        if frame is not None or update is not None:
            yield entry.metadata["timestamp"], frame, update
