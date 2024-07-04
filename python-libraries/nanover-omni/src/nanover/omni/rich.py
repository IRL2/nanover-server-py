from textual.app import App, ComposeResult
from textual.containers import Horizontal, Vertical, Grid
from textual.reactive import reactive
from textual.widgets import Button, Label, RadioSet

from nanover.omni import OmniRunner


class OmniTextualApp(App):
    CSS_PATH = "rich.tcss"
    BINDINGS = [("q", "quit", "Quit")]
    TITLE = "NanoVer Omni Server"

    simulation = reactive(None)

    def __init__(self, omni: OmniRunner):
        super().__init__()
        self.omni = omni

    def on_mount(self):
        self.set_interval(1 / 50, self.update)

    def update(self):
        running = self.omni.simulation is not None
        paused = self.omni.paused

        name = "None" if not running else self.omni.simulation.name
        self.query_one("#status", Label).update(f"Running: {name}")
        self.query_one("#controls").disabled = not running
        self.query_one("#play", Button).disabled = not paused
        self.query_one("#pause", Button).disabled = paused

    def on_radio_set_changed(self, event: RadioSet.Changed):
        self.omni.load(int(event.pressed.id[1:]))
        self.simulation = self.omni.simulation

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
            with Vertical():
                yield Label(id="status")
            with Grid(id="simulations"):
                for i, simulation in enumerate(self.omni.simulations):
                    yield Button(f'Load "{simulation.name}"', id=f"_{i}")
            with Horizontal(id="controls"):
                yield Button("Play", id="play")
                yield Button("Pause", id="pause")
                yield Button("Step", id="step")
                yield Button("Reset", id="reset")
            with Horizontal():
                yield Button("Quit", id="quit")
