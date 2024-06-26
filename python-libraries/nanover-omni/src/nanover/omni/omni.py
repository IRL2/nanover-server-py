from concurrent.futures import ThreadPoolExecutor, Future
from queue import Queue
from typing import Protocol, List, Optional

from nanover.app import NanoverImdApplication
from nanover.trajectory.frame_server import (
    LOAD_COMMAND_KEY,
    LIST_COMMAND_KEY,
    NEXT_COMMAND_KEY,
)


class Simulation(Protocol):
    name: str

    def run(self, app_server: NanoverImdApplication, cancel: Queue):
        pass


class OmniRunner:
    def __init__(self, app_server: NanoverImdApplication):
        self._app_server = app_server

        self.simulations: List[Simulation] = []
        self._simulation_index = 0
        self._simulation_counter = 0

        app_server.server.register_command(LOAD_COMMAND_KEY, self.load)
        app_server.server.register_command(NEXT_COMMAND_KEY, self.next)
        app_server.server.register_command(LIST_COMMAND_KEY, self.list)

        self._threads = ThreadPoolExecutor(max_workers=1)
        self._cancel: Optional[Queue] = None
        self._run_task: Optional[Future] = None

    def close(self):
        self.cancel_run()
        self.app_server.close()

    @property
    def app_server(self):
        return self._app_server

    @property
    def simulation(self):
        return self.simulations[self._simulation_index]

    def add_simulation(self, simulation: Simulation):
        self.simulations.append(simulation)

    def load(self, index: int):
        self.cancel_run()
        self._simulation_index = int(index) % len(self.simulations)
        self.run()

    def next(self):
        self.load(self._simulation_index + 1)

    def list(self):
        return {"simulations": [simulation.name for simulation in self.simulations]}

    def run(self):
        if self.is_running:
            raise RuntimeError("Already running on a thread!")

        self._cancel = Queue()
        self._run_task = self._threads.submit(
            self.simulation.run, self.app_server, self._cancel
        )

    def cancel_run(self):
        if self._cancel is None:
            return

        self._cancel.put("cancel")
        self._cancel = None
        self._run_task.result()
        self._run_task = None

    @property
    def is_running(self) -> bool:
        # ideally we'd just check _run_task.running(), but there can be a delay
        # between the task starting and hitting the running state.
        return self._run_task is not None and (
            self._run_task.cancelled() or self._run_task.done()
        )

    def __enter__(self):
        pass

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
