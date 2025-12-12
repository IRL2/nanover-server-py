import logging
from concurrent.futures import ThreadPoolExecutor, Future
from contextlib import suppress
from queue import Queue, Empty
from ssl import SSLContext
from typing import List, Set

from .imd_app import NanoverImdApplication
from .selection import RenderingSelection
from nanover.core import AppServer, basic_info_string, Simulation
from nanover.imd.imd_force import InvalidInteractionError
from nanover.trajectory import FrameData, keys
from nanover.utilities.change_buffers import DictionaryChange
from nanover.utilities.timing import VariableIntervalGenerator

CLEAR_PREFIXES = {"avatar.", "play-area.", "selection.", "scene", "interaction."}


class OmniRunner:
    """
    Provides a NanoVer server that supports switching between multiple simulations.
    """

    @classmethod
    def with_basic_server(
        cls,
        *simulations: Simulation,
        name: str | None = None,
        address: str | None = None,
        port: int | None = None,
        ssl: SSLContext | None = None,
    ):
        """
        Construct this using a basic NanoVer server and an optional list of initial simulations.

        :param simulations: List of starting simulations to make available
        :param name: Optional server name to broadcast
        :param address: Optional server address to use
        :param port: Optional server port to use
        :param ssl: Optional SSL context to use for the WebSocket server
        """
        app_server = NanoverImdApplication.basic_server(
            name=name,
            address=address,
            port=port,
            ssl=ssl,
        )

        omni = cls(app_server)
        for simulation in simulations:
            omni.add_simulation(simulation)
        return omni

    @classmethod
    def from_client(
        cls,
        client: NanoverImdApplication,
        *simulations: Simulation,
    ):
        """
        Use an existing client as the endpoint for sending simulated frames, registering play/pause etc commands,
        and the source of state for iMD interactions.
        """
        omni = cls(client)
        for simulation in simulations:
            omni.add_simulation(simulation)
        return omni

    def __init__(self, app_server: AppServer):
        self._app_server = app_server

        self.simulations: List[Simulation] = []
        self._simulation_index = 0
        self.simulation_selections: dict[Simulation, Set[RenderingSelection]] = {}

        app_server.register_command(keys.LOAD_COMMAND, self.load)
        app_server.register_command(keys.NEXT_COMMAND, self.next)
        app_server.register_command(keys.LIST_COMMAND, self.list)

        app_server.register_command(keys.RESET_COMMAND, self.reset)
        app_server.register_command(keys.PAUSE_COMMAND, self.pause)
        app_server.register_command(keys.PLAY_COMMAND, self.play)
        app_server.register_command(keys.STEP_COMMAND, self.step)

        self._threads = ThreadPoolExecutor(max_workers=1)
        self._runner: InternalRunner | None = None
        self._run_task: Future | None = None

        self.failed_simulations: Set[Simulation] = set()
        self.logging = logging.getLogger(__name__)

    def close(self):
        """
        Stop simulations and shut down server.
        """
        self._cancel_run()
        self.app_server.close()

    def print_basic_info(self):
        """
        Print out basic runner info to the terminal.
        """
        print(basic_info_string(self.app_server))

        list = "\n".join(
            f'{index}: "{simulation.name}"'
            for index, simulation in enumerate(self.simulations)
        )
        print(f"Available simulations:\n{list}")

    @property
    def app_server(self):
        """
        The NanoVer application server used by this runner.
        """
        return self._app_server

    @property
    def simulation(self):
        """
        The currently selected simulation.
        """
        try:
            return self.simulations[self._simulation_index]
        except IndexError:
            return None

    @property
    def is_paused(self):
        """
        :return: True if the current simulation is paused, False otherwise.
        """
        return self._runner.is_paused if self._runner is not None else None

    @property
    def runner(self):
        return self._runner

    def add_simulation(self, simulation: Simulation):
        """
        Add a simulation to list of available simulations to switch between.
        :param simulation: The simulation to add
        """
        self.simulations.append(simulation)

    def set_simulation_selections(
        self, simulation: Simulation, *selections: RenderingSelection
    ):
        existing_selections = self.simulation_selections.get(simulation) or set()
        existing_selections.update(selections)
        self.simulation_selections[simulation] = existing_selections

    def _clear_state(self):
        with self.app_server.lock_state() as state:
            removals = {
                key
                for key in state.keys()
                if any(key.startswith(prefix) for prefix in CLEAR_PREFIXES)
            }
        self.app_server.clear_locks()
        self.app_server.update_state(DictionaryChange(removals=removals))

    def _load_simulation_selections(self):
        with self.app_server.lock_state() as state:
            next_selections = {
                selection.selection_id: selection.to_dictionary()
                for selection in self.simulation_selections.get(self.simulation, [])
            }
            prev_selections = {
                key
                for key in state.keys()
                if key.startswith("selection.")
                if key not in next_selections
            }
            change = DictionaryChange(
                updates=next_selections,
                removals=prev_selections,
            )
        self.app_server.update_state(change)

    def load(self, index: int):
        """
        Switch to the simulation at a given index.
        :param index: Index of simulation to switch to
        :return:
        """
        self._cancel_run()
        self._simulation_index = int(index) % len(self.simulations)
        self._clear_state()
        self._load_simulation_selections()
        self._start_run()

    def next(self):
        """
        Switch to the next simulation in the list of available simulations, looping after the last.
        """
        self.load(self._simulation_index + 1)

    def list(self):
        """
        Get the list of available simulations.
        """
        return {"simulations": [simulation.name for simulation in self.simulations]}

    def reset(self):
        """
        Reset the currently-active simulation to its initial state.
        """
        assert self._runner is not None
        self._runner.signals.put("reset")

    def pause(self):
        """
        Pause the currently-active simulation.
        """
        assert self._runner is not None
        self._runner.signals.put("pause")

    def play(self):
        """
        Unpause the currently-active simulation.
        """
        assert self._runner is not None
        self._runner.signals.put("play")

    def step(self):
        """
        Step to the next frame in the currently-active simulation.
        """
        assert self._runner is not None
        self._runner.signals.put("step")

    def _start_run(self):
        if self._run_task is not None:
            raise RuntimeError("Already running on a thread!")

        self._runner = InternalRunner(self, self.simulation, self.app_server)
        self._run_task = self._threads.submit(self._runner.run)

    def _cancel_run(self):
        if self._runner is not None:
            self._runner.signals.put("cancel")
            self._runner = None

        if self._run_task is not None:
            with suppress(Exception):
                try:
                    self._run_task.result(timeout=5)
                except TimeoutError:
                    self.logging.warning("Timed out waiting for run task.")
            self._run_task = None

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()


class InternalRunner:
    def __init__(
        self,
        omni: OmniRunner,
        simulation: Simulation,
        app_server: AppServer,
    ):
        self.omni = omni
        self.simulation = simulation
        self.app_server = app_server

        self.signals: Queue[str] = Queue()
        self.cancelled = False
        self.is_paused = False

        self.variable_interval_generator = VariableIntervalGenerator(1 / 30)
        self.logger = logging.getLogger(simulation.name)

    @property
    def play_step_interval(self):
        return self.variable_interval_generator.interval

    @play_step_interval.setter
    def play_step_interval(self, interval: float):
        self.variable_interval_generator.interval = interval

    def run(self):
        def send_exception_frame(message):
            frame = FrameData()
            frame.simulation_exception = message
            self.app_server.frame_publisher.send_frame(frame)
            self.logger.exception(message)

        try:
            self.simulation.load()
            self.simulation.reset(self.app_server)
            self.omni.failed_simulations.discard(self.simulation)

            for dt in self.variable_interval_generator.yield_interval():
                try:
                    self.handle_signals()

                    if self.cancelled:
                        break
                    if not self.is_paused:
                        # for recording playback we want to know real time elapsed, for live simulations it is typically
                        # ignored and stepped one frame per invocation
                        self.simulation.advance_by_seconds(dt)
                except InvalidInteractionError:
                    send_exception_frame(
                        f"Invalid interaction, tried erasing all interactions during simulation `{self.simulation.name}`"
                    )
                    self.app_server.imd.clear_interactions()
                except Exception as e:
                    send_exception_frame(
                        f"{type(e).__name__} during simulation `{self.simulation.name}`"
                    )

                    self.is_paused = True
                    self.logger.warning("Simulation paused due to exception.")

        except Exception as e:
            self.omni.failed_simulations.add(self.simulation)
            self.logger.exception(
                f"{type(e)} loading simulation `{self.simulation.name}`"
            )
            self.app_server.frame_publisher.send_clear()
            raise

    def handle_signals(self):
        with suppress(Empty):
            while signal := self.signals.get_nowait():
                match signal:
                    case "pause":
                        self.is_paused = True
                    case "play":
                        self.is_paused = False
                    case "step":
                        self.is_paused = True
                        self.simulation.advance_by_one_step()
                    case "reset":
                        self.simulation.reset(self.app_server)
                    case "cancel":
                        self.cancelled = True
