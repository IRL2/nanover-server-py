from nanover.app import OmniRunner
from nanover.utilities.change_buffers import DictionaryChange
from nanover.websocket.record import record_from_runner, BackgroundRecordingContext
from nanover.imd.imd_state import ParticleInteraction


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
