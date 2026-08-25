import logging
from collections.abc import Iterable, Sequence
from functools import partial
from itertools import count
from typing import Any, TypedDict

import numpy as np
from ipywidgets import Output

from nanover.app import OmniRunner
from nanover.app.selection import (
    INTERACTION_METHOD_DEFAULT,
    KEY_PROPERTY_HIDE,
    KEY_PROPERTY_INTERACTION_METHOD,
    KEY_PROPERTY_RENDERER,
    KEY_PROPERTY_VELOCITY_RESET,
    KEY_SELECTED_PARTICLE_IDS,
    RENDERER_DEFAULT,
)
from nanover.core.app_server import StateService
from nanover.core.types import CommandHandler
from nanover.imd import ParticleInteraction
from nanover.imd.imd_state import interaction_to_dict
from nanover.trajectory import FrameData
from nanover.utilities.change_buffers import DictionaryChange
from nanover.utilities.transforms import Transform, matrix_from_state_transform
from nanover.websocket.client.app_client import NanoverImdClient
from nanover.websocket.record import BackgroundRecordingContext, record_from_runner

from .modes import Mode


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
        self.panels = PanelsUtility(runner.app_server)
        self.interactions = InteractionsUtility(runner.app_server)
        self.selections = SelectionsUtility(runner.app_server)
        self.transforms = TransformsUtility(runner.app_server)
        self.handles = TransformHandlesUtility(runner.app_server)
        self.modes = ModesManager(self)

    @property
    def scene_transform(self) -> Transform:
        return self.transforms.fetch_transform("simulation") or Transform.identity()

    def show_logging(self):
        output = Output()

        class Handler(logging.Handler):
            def emit(self, record: logging.LogRecord):
                with output:
                    print(self.format(record))

        logging.getLogger().addHandler(Handler())
        return output

    def notify_all(self, message: str):
        for command in self.runner.app_server.commands:
            if command.endswith("/notify"):
                self.runner.app_server.run_command(command, {"message": message})

    def get_shared_state_value(self, key: str, default=None):
        with self.runner.app_server.lock_state() as state:
            return state.get(key, default)

    def set_shared_state_value(self, key: str, value):
        self.runner.app_server.update_state(DictionaryChange(updates={key: value}))

    def define_command(
        self, key: str, *, handler: CommandHandler, icon="❓", label="unnamed command"
    ):
        self.runner.app_server.register_command(
            key, callback=handler, icon=icon, label=label
        )

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
        self.define_command(
            "user/recording/start",
            handler=self.start_recording,
            icon="⏺️",
            label="start recording",
        )
        self.define_command(
            "user/recording/stop",
            handler=self.stop_recording,
            icon="⏹️",
            label="stop recording",
        )
        self.define_command(
            "user/recording/checkpoint",
            handler=self.mark_checkpoint,
            icon="🚩",
            label="mark checkpoint",
        )

    def add_interaction_mode(self, mode: Mode, name: str, icon="👆"):
        # handle old case
        import inspect

        if inspect.isclass(mode):
            import warnings

            mode = mode()
            warnings.warn(
                "Interaction modes should be instances not classes",
                DeprecationWarning,
                stacklevel=2,
            )

        self.modes.add_mode(mode, name, icon)

    def use_transform_handles(self):
        from .transform_handles import use_transform_handles

        use_transform_handles(self)

    def intersect_transform_handles(self, point) -> dict | None:
        for key, handle in self.handles.all_prefixed_items():
            object_to_root = self.transforms.fetch_transform_root(handle["parent"])
            local_point = object_to_root.points_parent_to_local(point)
            center, radius = handle["sphere"]
            if np.linalg.norm(np.subtract(local_point, center)) < radius:
                return handle
        return None

    def make_ghost_from_mdanalysis(self, key: str, *, atoms):
        from nanover.jupyter.ghosts import GhostMoleculeData, GhostMoleculeObject

        return GhostMoleculeObject.from_ghost_data(
            key,
            ghost_data=GhostMoleculeData.from_atom_group(atoms),
            utilities=self,
        )

    def make_ghost_from_frame_data(
        self, key: str, *, frame_data: FrameData, atom_indices: Iterable[int]
    ):
        from nanover.jupyter.ghosts import GhostMoleculeData, GhostMoleculeObject

        ghost_data = GhostMoleculeData.from_frame_data(
            frame_data, atom_indices=atom_indices
        )

        return GhostMoleculeObject.from_ghost_data(
            key,
            ghost_data=ghost_data,
            utilities=self,
        )


class ModesManager:
    def __init__(self, utilities: NanoverJupyterUtilities):
        self._utilities = utilities
        self._active_mode = Mode()
        self._modes = {"normal": Mode()}

        prev_cursors: dict[str, dict] = {}

        def on_interaction_started(*, key: str, interaction: ParticleInteraction):
            self._active_mode.on_interaction_started(key=key, interaction=interaction)

        def on_interaction_stopped(*, key: str, interaction: ParticleInteraction):
            self._active_mode.on_interaction_stopped(key=key, interaction=interaction)

        def on_interaction_updated(*, key: str, interaction: ParticleInteraction):
            self._active_mode.on_interaction_updated(key=key, interaction=interaction)

        def on_cursor_updated(*, key: str, cursor: dict):
            prev_cursor = prev_cursors.get(key, {})
            prev_cursors[key] = cursor

            self._active_mode.on_cursor_updated(key=key, cursor=cursor)

            prev_held = set(prev_cursor.get("heldbuttons", []))
            next_held = set(cursor.get("heldbuttons", []))
            pressed = next_held - prev_held
            released = prev_held - next_held

            for button in pressed:
                self._active_mode.on_button_pressed(
                    key=key, cursor=cursor, button=button
                )
            for button in released:
                self._active_mode.on_button_released(
                    key=key, cursor=cursor, button=button
                )

        self._utilities.runner.app_server.imd.interaction_started.add_callback(
            on_interaction_started
        )
        self._utilities.runner.app_server.imd.interaction_stopped.add_callback(
            on_interaction_stopped
        )

        self._utilities.runner.app_server.imd.interaction_updated.add_callback(
            on_interaction_updated
        )

        def on_state_updated(*, access_token: str, change: DictionaryChange):
            for key, value in change.updates.items():
                if key.startswith("cursor"):
                    on_cursor_updated(key=key, cursor=value)

        self._utilities.runner.app_server.state_dictionary.content_updated.add_callback(
            on_state_updated
        )

    def enter_mode(self, name="normal", *, notify=True):
        if notify:
            self._utilities.notify_all(f"MODE {name}")
        self._active_mode.on_exit()
        self._active_mode = self._modes[name]
        self._active_mode.on_enter()

    def add_normal_mode(self):
        self.add_mode(Mode(), "normal")

    def add_mode(self, mode: Mode, name: str, icon="👆"):
        self._modes[name] = mode
        self._utilities.define_command(
            f"user/mode/{name}",
            handler=lambda: self.enter_mode(name),
            icon=icon,
            label=f"{name} mode",
        )

    def remove_mode(self, name: str):
        mode = self._modes.pop(name, None)

        if mode is not None:
            self._utilities.runner.app_server.unregister_command(f"user/mode/{name}")

    def clear_all_modes(self):
        for mode in set(self._modes.keys()):
            self.remove_mode(mode)

        self._modes.clear()
        self._modes["normal"] = Mode()
        self.enter_mode()


class StateKeysUtility:
    prefix = ""

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

    def all_prefixed(self):
        with self._state.lock_state() as state:
            return {
                key.removeprefix(self.prefix)
                for key in state
                if key.startswith(self.prefix)
            }

    def all_prefixed_items(self):
        with self._state.lock_state() as state:
            return {
                key.removeprefix(self.prefix): value
                for key, value in state.items()
                if key.startswith(self.prefix)
            }.items()

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

    def clear_all(self):
        self._buffer = DictionaryChange(removals=self.all_prefixed())
        self._keys = set()
        self.check_flush()


class SelectionsUtility(StateKeysUtility):
    prefix = "selection."

    def update_selection(
        self,
        key: str,
        *,
        particle_ids: list[int] | None = None,
        renderer=RENDERER_DEFAULT,
        interaction_method=INTERACTION_METHOD_DEFAULT,
        velocity_reset=False,
        hide=False,
    ):
        self.update_object(
            f"{self.prefix}{key}",
            {
                "id": f"{self.prefix}{key}",
                "selected": {
                    KEY_SELECTED_PARTICLE_IDS: particle_ids
                    if particle_ids is not None
                    else [],
                },
                "properties": {
                    KEY_PROPERTY_RENDERER: renderer,
                    KEY_PROPERTY_INTERACTION_METHOD: interaction_method,
                    KEY_PROPERTY_VELOCITY_RESET: velocity_reset,
                    KEY_PROPERTY_HIDE: hide,
                },
            },
        )

    def remove_selection(self, key: str):
        self.remove_object(f"{self.prefix}{key}")


class PanelsUtility(StateKeysUtility):
    prefix = "panel."

    @staticmethod
    def header(
        label="header",
    ):
        return {
            "type": "header",
            "label": label,
        }

    @staticmethod
    def button(
        label="button",
        command="test/hello",
        arguments=None,
    ):
        if arguments is None:
            arguments = {}
        return {
            "type": "button",
            "label": label,
            "command": command,
            "arguments": arguments,
        }

    @staticmethod
    def slider(
        label="slider",
        variable="variable.dummy",
        range=(0.0, 1.0),
        integer=False,
        step=None,
    ):
        return {
            "type": "slider",
            "label": label,
            "variable": variable,
            "range": range,
            "integer": integer,
            "step": step,
        }

    def update_panel(
        self,
        key: str,
        *elements: Any,
        position=(0.0, 0.0, 0.0),
        label="Unnamed panel",
        **kwargs,
    ):
        self.update_object(
            f"{self.prefix}{key}",
            {
                "position": position,
                "label": label,
                "elements": elements,
                **kwargs,
            },
        )

    def remove_panel(self, key: str):
        self.remove_object(f"{self.prefix}{key}")


class InteractionsUtility(StateKeysUtility):
    prefix = "interaction."

    def update_interaction(self, key: str, interaction: ParticleInteraction):
        self.update_object(f"{self.prefix}{key}", interaction_to_dict(interaction))

    def remove_interaction(self, key: str):
        self.remove_object(f"{self.prefix}{key}")


class StateTransformEntry(TypedDict):
    transform: Sequence[float]
    parent: str | None


class TransformsUtility(StateKeysUtility):
    prefix = "transform."

    def update_transform(
        self,
        key: str,
        *,
        transform: Transform,
        parent="simulation",
    ):
        self.update_object(
            f"{self.prefix}{key}",
            {
                "transform": transform.to_state_transform(),
                "parent": parent,
            },
        )

    def remove_transform(
        self,
        key: str,
    ):
        self.remove_object(f"{self.prefix}{key}")

    def get_parent(self, key: str, *, default=None):
        entry = self.fetch_transform_entry(key)
        return default if entry is None else entry.get("parent", default)

    def fetch_transform_entry(self, key: str) -> StateTransformEntry | None:
        with self._state.lock_state() as state:
            return state.get(f"{self.prefix}{key}", None)

    def fetch_transform(self, key: str, *, default: Transform | None = None):
        entry = self.fetch_transform_entry(key)
        return (
            Transform.from_state_transform(entry["transform"])
            if entry is not None
            else default
        )

    def fetch_transform_root(self, key: str):
        matrix = np.identity(4)

        with self._state.lock_state() as state:
            while key is not None:
                entry: StateTransformEntry | None = state.get(
                    f"{self.prefix}{key}", None
                )

                if not entry:
                    break

                matrix = matrix_from_state_transform(entry["transform"]) @ matrix
                key = entry["parent"]

        return Transform.from_local_to_parent_matrix(matrix)


class TransformHandlesUtility(StateKeysUtility):
    prefix = "handle."

    def update_handle(
        self,
        key: str,
        *,
        parent: str,
        sphere=((0, 0, 0), 0.25),
        translate=True,
        rotate=True,
        scale=False,
        **kwargs,
    ):
        self.update_object(
            f"{self.prefix}{key}",
            dict(
                parent=parent,
                sphere=sphere,
                translate=translate,
                rotate=rotate,
                scale=scale,
                **kwargs,
            ),
        )

    def remove_handle(self, key: str):
        self.remove_object(f"{self.prefix}{key}")


class SceneObjectsUtility(StateKeysUtility):
    prefix = "object."

    def update_shape(
        self,
        key: str,
        *,
        shape="sphere",
        position=(0.0, 0.0, 0.0),
        color=(1.0, 1.0, 1.0, 1.0),
        size=0.1,
        parent="simulation",
        **kwargs,
    ):
        self.update_object(
            f"object.shape.{key}",
            {
                "shape": shape,
                "position": position,
                "color": color,
                "size": size,
                "parent": parent,
                **kwargs,
            },
        )

    def update_line(
        self,
        key: str,
        *,
        positions=((0.0, 0.0, 0.0), (1.0, 1.0, 1.0)),
        colors=None,
        color=(1.0, 1.0, 1.0, 1.0),
        size=0.05,
        parent="simulation",
        **kwargs,
    ):
        self.update_object(
            f"object.line.{key}",
            {
                "positions": positions,
                "colors": colors,
                "color": color,
                "size": size,
                "parent": parent,
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
        parent="simulation",
        **kwargs,
    ):
        self.update_object(
            f"object.label.{key}",
            {
                "text": text,
                "position": position,
                "color": color,
                "size": size,
                "parent": parent,
                **kwargs,
            },
        )

    def update_sprite(
        self,
        key: str,
        *,
        texture: str,
        position=(0.0, 0.0, 0.0),
        color=(1.0, 1.0, 1.0, 1.0),
        size=1.0,
        parent="simulation",
        **kwargs,
    ):
        self.update_object(
            f"object.sprite.{key}",
            dict(
                texture=texture,
                position=position,
                color=color,
                size=size,
                parent=parent,
                **kwargs,
            ),
        )

    def remove_shape(self, key: str):
        self.remove_object(f"object.shape.{key}")

    def remove_line(self, key: str):
        self.remove_object(f"object.line.{key}")

    def remove_label(self, key: str):
        self.remove_object(f"object.label.{key}")


def make_id_generator(prefix=""):
    return partial(next, (f"{prefix}{i}" for i in count()))
