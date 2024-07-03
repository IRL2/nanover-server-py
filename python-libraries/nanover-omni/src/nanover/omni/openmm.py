from pathlib import Path
from typing import Optional, Any

from openmm.app import Simulation

from nanover.app import NanoverImdApplication
from nanover.openmm import serializer
from nanover.openmm.imd import create_imd_force, NanoverImdReporter


class OpenMMSimulation:
    def __init__(self, path: str):
        self.name = Path(path).stem

        self.xml_path = path
        self.app_server: Optional[NanoverImdApplication] = None

        self.frame_interval = 5
        self.force_interval = 5

        self.imd_force = create_imd_force()
        self.simulation: Optional[Simulation] = None
        self.checkpoint: Optional[Any] = None
        self.reporter: Optional[NanoverImdReporter] = None

    def load(self):
        platform = None

        with open(str(self.xml_path)) as infile:
            self.imd_force = create_imd_force()
            self.simulation = serializer.deserialize_simulation(
                infile.read(), imd_force=self.imd_force, platform_name=platform
            )

        self.checkpoint = self.simulation.context.createCheckpoint()

    def unload(self):
        self.simulation = None
        self.checkpoint = None

    def reset(self, app_server: NanoverImdApplication):
        assert self.simulation is not None and self.checkpoint is not None

        self.app_server = app_server
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
        )
        self.simulation.reporters.append(self.reporter)

    def advance_to_next_report(self):
        assert self.simulation is not None
        self.simulation.step(self.frame_interval)

    def advance_by_seconds(self, dt: float):
        self.advance_to_next_report()

    def advance_by_one_step(self):
        self.advance_to_next_report()
