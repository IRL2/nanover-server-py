from ..trajectory.frame_server import FrameServer
from .runner import Runner
from .narupareporter import NarupaReporter


class Server(Runner):
    def __init__(self, simulation, *, address: str, port: int, reportInterval: int=1):
        super().__init__(simulation)
        self._frame_server = FrameServer(address=address, port=port)
        self._frame_reporter = NarupaReporter(
            reportInterval=reportInterval,
            frameServer=self._frame_server,
        )
        self.make_report()

    def make_report(self):
        if not self.reporting:
            self.simulation.reporters.append(self._frame_reporter)

    def make_not_report(self):
        if self.reporting:
            self.simulation.reporters.remove(self._frame_reporter)

    @property
    def reporting(self):
        return self._frame_reporter in self.simulation.reporters

    @reporting.setter
    def reporting(self, value):
        if value:
            self.make_report()
        else:
            self.make_quiet()
