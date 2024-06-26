from pathlib import Path
from queue import Queue
from typing import Optional

import openmm

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

        self.frame_index = 0

        self.frame_interval = 5

    def run(self, app_server: NanoverImdApplication, cancel: Queue):
        platform = None

        with open(str(self.xml_path)) as infile:
            imd_force = create_imd_force()
            simulation = serializer.deserialize_simulation(
                infile.read(), imd_force=imd_force, platform_name=platform
            )

        simulation_counter = 0
        reporter = NanoverImdReporter(
            frame_interval=5,
            force_interval=5,
            imd_force=imd_force,
            imd_state=app_server.imd,
            frame_publisher=app_server.frame_publisher,
            simulation_counter=simulation_counter,
        )
        simulation.reporters.append(reporter)

        steps: Optional[int] = None
        remaining_steps = steps if steps is not None else float("inf")
        for _ in self._variable_interval_generator.yield_interval():
            if not cancel.empty() or remaining_steps <= 0:
                break
            steps_for_this_iteration = min(self.frame_interval, remaining_steps)
            try:
                simulation.step(steps_for_this_iteration)
            except (ValueError, openmm.OpenMMException):
                # We want to stop running if the simulation exploded in a way
                # that prevents OpenMM to run. Otherwise, we will be a state
                # where OpenMM raises an exception which would make the runner
                # unusable. The OpenMMException is typically raised by OpenMM
                # itself when something is NaN; the ValueError is typically
                # raised by the StateReporter when the energy is NaN.
                break
            remaining_steps -= steps_for_this_iteration
