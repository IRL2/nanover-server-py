from nanover.recording import NanoverRecordingReader
from nanover.utilities.change_buffers import DictionaryChange


def iter_checkpoints(path, checkpoint_key="mark.checkpoint"):
    with NanoverRecordingReader.from_path(path) as reader:
        for event in reader.iter_max():
            if event.next_state_event is None:
                continue

            change = DictionaryChange.from_dict(event.next_state_event.message)

            if checkpoint_key in change.updates:
                yield event
