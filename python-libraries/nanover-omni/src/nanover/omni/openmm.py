from pathlib import Path
from queue import Queue
from typing import Optional, Any

from openmm.app import Simulation

from nanover.app import NanoverImdApplication
from nanover.openmm import serializer
from nanover.openmm.imd import create_imd_force, NanoverImdReporter
from nanover.utilities.timing import VariableIntervalGenerator


class OpenMMSimulation:
    def __init__(self, path: str):
        self.name = Path(path).stem

        self.xml_path = path
        self.app_server: Optional[NanoverImdApplication] = None

        self._variable_interval_generator = VariableIntervalGenerator(1 / 30)

        self.frame_interval = 5
        self.force_interval = 5

        self.paused = False

        self.imd_force = create_imd_force()
        self.simulation: Optional[Simulation] = None
        self.checkpoint: Optional[Any] = None
        self.reporter: Optional[NanoverImdReporter] = None

    def load(self):
        platform = None

        with open(str(self.xml_path)) as infile:
            self.simulation = serializer.deserialize_simulation(
                infile.read(), imd_force=self.imd_force, platform_name=platform
            )

        self.checkpoint = self.simulation.context.createCheckpoint()

    def unload(self):
        self.simulation = None
        self.checkpoint = None

    def reset(self, simulation_counter=0):
        assert self.simulation is not None and self.checkpoint is not None

        self.simulation.context.loadCheckpoint(self.checkpoint)

        try:
            self.simulation.reporters.remove(self.reporter)
        except ValueError:
            pass

        self.reporter = NanoverImdReporter(
            frame_interval=self.frame_interval,
            force_interval=self.force_interval,
            imd_force=self.imd_force,
            imd_state=self.app_server.imd,
            frame_publisher=self.app_server.frame_publisher,
            simulation_counter=simulation_counter,
        )
        self.simulation.reporters.append(self.reporter)

    def run(self, app_server: NanoverImdApplication, cancel: Queue):
        self.app_server = app_server

        self.load()
        self.reset()

        for _ in self._variable_interval_generator.yield_interval():
            if not cancel.empty():
                break
            if not self.paused:
                self.advance_to_next_report()

    def advance_to_next_report(self):
        assert self.simulation is not None
        self.simulation.step(self.frame_interval)
