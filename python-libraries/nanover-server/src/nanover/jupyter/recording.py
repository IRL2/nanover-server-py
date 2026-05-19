from nanover.app import OmniRunner
from nanover.websocket.record import record_from_runner, BackgroundRecordingContext


class RecordingCommands:
    @classmethod
    def from_runner(cls, runner: OmniRunner):
        return cls(runner.app_server)

    def __init__(self, runner: OmniRunner):
        self.path_template = "RECORDING-{count}-{simulation_name}.nanover.zip"

        self._runner = runner
        runner.app_server.register_command("user/recording/start", self.start_recording)
        runner.app_server.register_command("user/recording/stop", self.stop_recording)

        self._recordings = []
        self._recorder: BackgroundRecordingContext | None = None

    def start_recording(self):
        self.stop_recording()
        out_path = self.path_template.format(
            count=len(self._recordings),
            simulation_name=self._runner.simulation.name
        )
        self._recorder = record_from_runner(self._runner, out_path)
        self._recordings.append(out_path)

    def stop_recording(self):
        if self._recorder is not None:
            self._recorder.close()
            self._recorder = None
