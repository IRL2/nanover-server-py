from concurrent.futures import ThreadPoolExecutor, Future
from contextlib import suppress
from queue import Queue, Empty
from typing import Protocol, List, Optional

from nanover.app import NanoverImdApplication
from nanover.trajectory.frame_server import (
    LOAD_COMMAND_KEY,
    LIST_COMMAND_KEY,
    NEXT_COMMAND_KEY,
    RESET_COMMAND_KEY,
    PAUSE_COMMAND_KEY,
    PLAY_COMMAND_KEY,
    STEP_COMMAND_KEY,
)
from nanover.utilities.timing import VariableIntervalGenerator


class Simulation(Protocol):
    name: str

    def load(self): ...
    def reset(self, app_server: NanoverImdApplication): ...
    def advance_by_one_step(self): ...
    def advance_by_seconds(self, dt: float): ...


class OmniRunner:
    @classmethod
    def with_basic_server(
        cls,
        name: Optional[str] = None,
        address: Optional[str] = None,
        port: Optional[int] = None,
    ):
        app_server = NanoverImdApplication.basic_server(name, address, port)
        return cls(app_server)

    def __init__(self, app_server: NanoverImdApplication):
        self._app_server = app_server

        self.simulations: List[Simulation] = []
        self._simulation_index = 0
        self._simulation_counter = -1

        app_server.server.register_command(LOAD_COMMAND_KEY, self.load)
        app_server.server.register_command(NEXT_COMMAND_KEY, self.next)
        app_server.server.register_command(LIST_COMMAND_KEY, self.list)

        app_server.server.register_command(RESET_COMMAND_KEY, self.reset)
        app_server.server.register_command(PAUSE_COMMAND_KEY, self.pause)
        app_server.server.register_command(PLAY_COMMAND_KEY, self.play)
        app_server.server.register_command(STEP_COMMAND_KEY, self.step)

        self._threads = ThreadPoolExecutor(max_workers=1)
        self._runner: Optional[InternalRunner] = None
        self._run_task: Optional[Future] = None

    def close(self):
        self.app_server.close()
        self._cancel_run()

    @property
    def app_server(self):
        return self._app_server

    @property
    def simulation(self):
        return self.simulations[self._simulation_index]

    def add_simulation(self, simulation: Simulation):
        self.simulations.append(simulation)

    def load(self, index: int):
        self._cancel_run()
        self._simulation_index = int(index) % len(self.simulations)
        self._simulation_counter += 1
        self._start_run()

    def next(self):
        self.load(self._simulation_index + 1)

    def list(self):
        return {"simulations": [simulation.name for simulation in self.simulations]}

    def reset(self):
        assert self._runner is not None
        self._simulation_counter += 1
        self._runner.signals.put("reset")

    def pause(self):
        assert self._runner is not None
        self._runner.signals.put("pause")

    def play(self):
        assert self._runner is not None
        self._runner.signals.put("play")

    def step(self):
        assert self._runner is not None
        self._runner.signals.put("step")

    def _start_run(self):
        if self._run_task is not None:
            raise RuntimeError("Already running on a thread!")

        self._runner = InternalRunner(self.simulation, self.app_server)
        self._run_task = self._threads.submit(self._runner.run)

    def _cancel_run(self):
        if self._runner is not None:
            self._runner.signals.put("cancel")
            self._runner = None

        if self._run_task is not None:
            self._run_task.result()
            self._run_task = None

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()


class InternalRunner:
    def __init__(self, simulation: Simulation, app_server: NanoverImdApplication):
        self.simulation = simulation
        self.app_server = app_server

        self.signals: Queue[str] = Queue()
        self.cancelled = False
        self.paused = False

        self.variable_interval_generator = VariableIntervalGenerator(1 / 30)

    def run(self):
        self.simulation.load()
        self.simulation.reset(self.app_server)

        for dt in self.variable_interval_generator.yield_interval():
            self.handle_signals()

            if self.cancelled:
                break
            if not self.paused:
                self.simulation.advance_by_seconds(dt)

    def handle_signals(self):
        with suppress(Empty):
            while signal := self.signals.get_nowait():
                match signal:
                    case "pause":
                        self.paused = True
                    case "play":
                        self.paused = False
                    case "step":
                        self.paused = True
                        self.simulation.advance_by_one_step()
                    case "reset":
                        self.simulation.reset(self.app_server)
                    case "cancel":
                        self.cancelled = True
