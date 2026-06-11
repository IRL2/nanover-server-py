from typing import Any

from nanover.app import OmniRunner
from nanover.core.app_server import StateService
from nanover.utilities.change_buffers import DictionaryChange
from nanover.websocket.client.app_client import NanoverImdClient
from nanover.websocket.record import record_from_runner, BackgroundRecordingContext
from nanover.imd.imd_state import (
    ParticleInteraction,
    interaction_to_dict,
    INTERACTION_PREFIX,
)


class Mode:
    def on_interaction_started(self, *, key: str, interaction: ParticleInteraction):
        pass

    def on_interaction_stopped(self, *, key: str, interaction: ParticleInteraction):
        pass


class NanoverJupyterUtilities:
    _recording_path: str | None = None
    _recording_count = 0
    _recorder: BackgroundRecordingContext | None = None
    _next_checkpoint_index = 0
    _active_mode = Mode()

    @classmethod
    def from_runner(cls, runner: OmniRunner):
        return cls(runner)

    def __init__(self, runner: OmniRunner):
        self.runner = runner
        self.objects = SceneObjectsUtility(runner.app_server)
        self.interactions = InteractionsUtility(runner.app_server)

    def notify_all(self, message: str):
        for command in self.runner.app_server.commands:
            if command.endswith("/notify"):
                self.runner.app_server.run_command(command, dict(message=message))

    def start_recording(self):
        self._recording_path = f"RECORDING-{self._recording_count}-{self.runner.simulation.name}.nanover.zip"
        self._recorder = record_from_runner(self.runner, self._recording_path)
        self._recording_count += 1
        self.notify_all(f"STARTED RECORDING to {self._recording_path}")

    def stop_recording(self):
        if self._recorder is not None:
            self._recorder.close()
            self.notify_all(f"FINISHED RECORDING to {self._recording_path}")

    def mark_checkpoint(self):
        self.runner.app_server.update_state(
            DictionaryChange(updates={"mark.checkpoint": self._next_checkpoint_index})
        )
        self.notify_all(f"MARKED CHECKPOINT {self._next_checkpoint_index}")
        self._next_checkpoint_index += 1

    def use_recording_commands(self):
        self.runner.app_server.register_command(
            "user/recording/start", self.start_recording
        )
        self.runner.app_server.register_command(
            "user/recording/stop", self.stop_recording
        )
        self.runner.app_server.register_command(
            "user/recording/checkpoint", self.mark_checkpoint
        )

    def use_interaction_modes(self):
        def on_interaction_started(*, key: str, interaction: ParticleInteraction):
            self._active_mode.on_interaction_started(key=key, interaction=interaction)

        def on_interaction_stopped(*, key: str, interaction: ParticleInteraction):
            self._active_mode.on_interaction_stopped(key=key, interaction=interaction)

        self.runner.app_server.imd.interaction_started.add_callback(
            on_interaction_started
        )
        self.runner.app_server.imd.interaction_stopped.add_callback(
            on_interaction_stopped
        )

        self._active_mode = Mode()
        self.add_interaction_mode(Mode, "normal")

    def add_interaction_mode[T: type[Mode]](self, mode: T, name: str):
        def enter():
            self._active_mode = mode()
            self.notify_all(f"INTERACTION MODE {name}")

        self.runner.app_server.register_command(f"user/interaction/{name}", enter)


class StateKeysUtility:
    @classmethod
    def from_runner(cls, runner: OmniRunner):
        return cls(runner.app_server)

    @classmethod
    def from_client(cls, client: NanoverImdClient):
        return cls(client)

    def __init__(self, state: StateService):
        self._state = state
        self._buffer = DictionaryChange()
        self._keys: set[str] = set()
        self._depth = 0

    def __enter__(self):
        self._depth += 1
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self._depth -= 1
        self.check_flush()

    def update_object(self, key: str, value: dict[str, Any]):
        self._buffer.updates[key] = value
        self.check_flush()

    def remove_object(self, key: str):
        self._buffer.updates.pop(key, None)
        self._buffer.removals = {key, *self._buffer.removals}
        self.check_flush()

    def check_flush(self):
        if self._depth == 0:
            self.flush()

    def flush(self):
        self._keys |= self._buffer.updates.keys()
        self._state.update_state(self._buffer)
        self._buffer = DictionaryChange()

    def clear(self):
        self._buffer = DictionaryChange(removals=self._keys)
        self._keys = set()
        self.check_flush()


class InteractionsUtility(StateKeysUtility):
    def clear_all(self):
        keys = {
            key
            for key in self._state.state_dictionary.copy_content()
            if key.startswith("interaction.")
        }
        self._buffer = DictionaryChange(removals=keys)
        self._keys = set()
        self.check_flush()

    def update_interaction(self, key: str, interaction: ParticleInteraction):
        self.update_object(
            f"{INTERACTION_PREFIX}{key}", interaction_to_dict(interaction)
        )

    def remove_interaction(self, key: str):
        self.remove_object(f"{INTERACTION_PREFIX}{key}")


class SceneObjectsUtility(StateKeysUtility):
    def clear_all(self):
        keys = {
            key
            for key in self._state.state_dictionary.copy_content()
            if key.startswith("object.")
        }
        self._buffer = DictionaryChange(removals=keys)
        self._keys = set()
        self.check_flush()

    def update_shape(
        self,
        key: str,
        *,
        shape="sphere",
        position=(0.0, 0.0, 0.0),
        color=(1.0, 1.0, 1.0, 1.0),
        size=0.1,
        **kwargs,
    ):
        self.update_object(
            f"object.shape.{key}",
            {
                "shape": shape,
                "position": position,
                "color": color,
                "size": size,
                **kwargs,
            },
        )

    def update_line(
        self,
        key: str,
        *,
        positions=((0.0, 0.0, 0.0), (1.0, 1.0, 1.0)),
        color=(1.0, 1.0, 1.0, 1.0),
        size=0.05,
        **kwargs,
    ):
        self.update_object(
            f"object.line.{key}",
            {
                "positions": positions,
                "color": color,
                "size": size,
                **kwargs,
            },
        )

    def update_label(
        self,
        key: str,
        *,
        text="label",
        position=(0.0, 0.0, 0.0),
        color=(1.0, 1.0, 1.0, 1.0),
        size=0.05,
        **kwargs,
    ):
        self.update_object(
            f"object.label.{key}",
            {
                "text": text,
                "position": position,
                "color": color,
                "size": size,
                **kwargs,
            },
        )

    def remove_shape(self, key: str):
        self.remove_object(f"object.shape.{key}")

    def remove_line(self, key: str):
        self.remove_object(f"object.line.{key}")

    def remove_label(self, key: str):
        self.remove_object(f"object.label.{key}")
