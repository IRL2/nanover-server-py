from textual.app import App, ComposeResult
from textual.containers import Horizontal, Vertical, VerticalScroll
from textual.widgets import Button, Label

from nanover.omni import OmniRunner


class OmniTextualApp(App):
    CSS_PATH = "rich.tcss"
    BINDINGS = [("q", "quit", "Quit")]
    TITLE = "NanoVer Omni Server"

    def __init__(self, omni: OmniRunner):
        super().__init__()
        self.omni = omni

    def on_mount(self):
        self.set_interval(1 / 2, self.update)

    def update(self):
        running = self.omni.simulation is not None
        paused = self.omni.is_paused
        failed = self.omni.simulation in self.omni.failed_simulations

        name = "None" if not running else self.omni.simulation.name
        self.query_one("#status Label", Label).update(
            f"Running: {name}" if not failed else "Failed."
        )
        self.query_one("#play", Button).disabled = not paused
        self.query_one("#pause", Button).disabled = paused

        for button in self.query("#controls Button"):
            button.disabled = not running
        self.query_one("#quit").disabled = False

        for i, simulation in enumerate(self.omni.simulations):
            failed = simulation in self.omni.failed_simulations
            self.query_one(f"#_{i}").classes = "failed" if failed else ""

    async def on_button_pressed(self, event: Button.Pressed):
        match event.button.id:
            case "play":
                self.omni.play()
            case "pause":
                self.omni.pause()
            case "step":
                self.omni.step()
            case "reset":
                self.omni.reset()
            case "quit":
                await self.run_action("quit")
            case button_id:
                self.omni.load(int(button_id[1:]))

    def compose(self) -> ComposeResult:
        with Vertical(id="container"):
            with Vertical(id="status"):
                yield Label()
            with VerticalScroll(id="simulations"):
                for i, simulation in enumerate(self.omni.simulations):
                    failed = simulation in self.omni.failed_simulations
                    yield Button(
                        f"{simulation.name}",
                        id=f"_{i}",
                        classes="failed" if failed else None,
                    )
            with Horizontal(id="controls"):
                with Horizontal():
                    yield Button("Play", id="play")
                    yield Button("Pause", id="pause")
                    yield Button("Step", id="step")
                    yield Button("Reset", id="reset")
                    yield Button("Quit", id="quit")
