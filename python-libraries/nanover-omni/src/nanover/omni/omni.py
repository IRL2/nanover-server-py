from typing import Protocol, List

from nanover.app import NanoverImdApplication
from nanover.trajectory.frame_server import (
    LOAD_COMMAND_KEY,
    LIST_COMMAND_KEY,
    NEXT_COMMAND_KEY,
)


class Simulation(Protocol):
    name: str

    def load(self, app_server: NanoverImdApplication):
        pass

    def run(self):
        pass

    def cancel_run(self):
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

    def close(self):
        self.simulation.cancel_run()
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
        self.simulation.cancel_run()
        self._simulation_index = int(index) % len(self.simulations)
        self.simulation.load(self.app_server)
        self.simulation.run()

    def next(self):
        self.load(self._simulation_index + 1)

    def list(self):
        return {"simulations": [simulation.name for simulation in self.simulations]}

    def __enter__(self):
        pass

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
