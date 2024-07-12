import logging
from os import PathLike
from pathlib import Path
from typing import Optional, Any

from openmm.app import Simulation, StateDataReporter

from nanover.app import NanoverImdApplication
from nanover.openmm import serializer
from nanover.openmm.imd import (
    create_imd_force,
    NanoverImdReporter,
    get_imd_forces_from_system,
)


class OpenMMSimulation:
    @classmethod
    def from_simulation(cls, simulation: Simulation, *, name: Optional[str] = None):
        sim = cls(name)

        sim.simulation = simulation
        potential_imd_forces = get_imd_forces_from_system(sim.simulation.system)
        if not potential_imd_forces:
            raise ValueError(
                "The simulation must include an appropriate force for imd."
            )
        if len(potential_imd_forces) > 1:
            logging.warning(
                f"More than one force could be used as imd force "
                f"({len(potential_imd_forces)}); taking the last one."
            )
        # In case there is more than one compatible force we take the last one.
        # The forces are in the order they have been added, so we take the last
        # one that have been added. This is the most likely to have been added
        # for the purpose of this runner, the other ones are likely leftovers
        # or forces created for another purpose.
        sim.imd_force = potential_imd_forces[-1]

        sim.checkpoint = sim.simulation.context.createCheckpoint()

        return sim

    @classmethod
    def from_xml_path(cls, path: PathLike[str], *, name: Optional[str] = None):
        sim = cls(name or Path(path).stem)
        sim.xml_path = path
        return sim

    def __init__(self, name: Optional[str] = None):
        self.name = name or "Unnamed OpenMM Simulation"

        self.xml_path: Optional[PathLike[str]] = None
        self.app_server: Optional[NanoverImdApplication] = None

        self.frame_interval = 5
        self.force_interval = 5
        self.platform: Optional[str] = None

        self.imd_force = create_imd_force()
        self.simulation: Optional[Simulation] = None
        self.checkpoint: Optional[Any] = None
        self.reporter: Optional[NanoverImdReporter] = None
        self.verbose_reporter: Optional[StateDataReporter] = None

    def load(self):
        if self.xml_path is None or self.simulation is not None:
            return

        with open(self.xml_path) as infile:
            self.imd_force = create_imd_force()
            self.simulation = serializer.deserialize_simulation(
                infile, imd_force=self.imd_force, platform_name=self.platform
            )

        self.checkpoint = self.simulation.context.createCheckpoint()

    def reset(self, app_server: NanoverImdApplication):
        assert self.simulation is not None and self.checkpoint is not None

        self.app_server = app_server
        self.simulation.context.loadCheckpoint(self.checkpoint)

        try:
            self.simulation.reporters.remove(self.reporter)
            if self.verbose_reporter is not None:
                self.simulation.reporters.remove(self.verbose_reporter)
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
        if self.verbose_reporter is not None:
            self.simulation.reporters.append(self.verbose_reporter)

    def advance_to_next_report(self):
        assert self.simulation is not None
        self.simulation.step(self.frame_interval)

    def advance_by_seconds(self, dt: float):
        self.advance_to_next_report()

    def advance_by_one_step(self):
        self.advance_to_next_report()
