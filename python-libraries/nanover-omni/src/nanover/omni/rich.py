from textual.app import App, ComposeResult
from textual.containers import Horizontal, Vertical, VerticalScroll
from textual.css.query import NoMatches
from textual.reactive import reactive
from textual.screen import Screen
from textual.widgets import Button, Label, Static

from nanover.omni import OmniRunner


class BSOD(Screen):
    BINDINGS = [("escape", "app.pop_screen", "Pop screen")]

    def compose(self) -> ComposeResult:
        yield Static(" Windows ", id="title")
        yield Static(ERROR_TEXT)
        yield Static("Press any key to continue [blink]_[/]", id="any-key")


class OmniTextualApp(App):
    CSS_PATH = "rich.tcss"
    BINDINGS = [("q", "quit", "Quit"), ("escape", "push_screen('bsod')", "BSOD")]
    TITLE = "NanoVer Omni Server"

    SCREENS = {"bsod": BSOD()}

    def __init__(self, omni: OmniRunner):
        super().__init__()
        self.omni = omni

    def on_mount(self):
        self.set_interval(1 / 50, self.update)

    def update(self):
        if not self.screen.is_current:
            return

        running = self.omni.simulation is not None
        paused = self.omni.paused

        # TODO: don't know how to properly avoid calling this when elements absent...
        try:
            name = "None" if not running else self.omni.simulation.name
            self.query_one("#status Label", Label).update(f"Running: {name}")
            self.query_one("#play", Button).disabled = not paused
            self.query_one("#pause", Button).disabled = paused

            for button in self.query("#controls Button"):
                button.disabled = not running
            self.query_one("#quit").disabled = False
        except NoMatches:
            ...

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
                    yield Button(f"{simulation.name}", id=f"_{i}")
            with Horizontal(id="controls"):
                with Horizontal():
                    yield Button("Play", id="play")
                    yield Button("Pause", id="pause")
                    yield Button("Step", id="step")
                    yield Button("Reset", id="reset")
                    yield Button("Quit", id="quit")


ERROR_TEXT = """
An error has occurred. To continue:

Press Enter to return to Windows, or

Press CTRL+ALT+DEL to restart your computer. If you do this,
you will lose any unsaved information in all open applications.

Error: 0E : 016F : BFF9B3D4
"""

