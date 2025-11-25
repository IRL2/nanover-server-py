from typing import Any

from ipywidgets import interact
import ipywidgets as widgets
from IPython.display import display
from nanover.core import AppServerMinimal

from nanover.omni import OmniRunner
from nanover.trajectory import keys
from nanover.utilities.change_buffers import DictionaryChange


def interact_button(on_click, *, description: str):
    button = make_button(on_click, description=description)
    display(button)
    return button


def make_button(on_click, *, description: str):
    button = widgets.Button(description=description)
    button.on_click(on_click)
    return button


def make_toggle(setter, *, description: str):
    toggle = widgets.ToggleButton(description=description)
    toggle.observe(lambda change: setter(change["new"]), names="value")
    return toggle


def show_runner_controls(imd_runner: OmniRunner):
    def close_server(_):
        imd_runner.close()

    interact_button(close_server, description="Close Server"),

    show_app_server_controls(imd_runner.app_server)


def show_app_server_controls(app_server: AppServerMinimal):
    copied_box: Any = None

    def get_shared_value(key: str, default: Any = None):
        return app_server.copy_state().get(key, default)

    def set_shared_value(key: str, value: Any):
        updates = {key: value}
        change = DictionaryChange(updates=updates)
        app_server.update_state(change)

    def make_state_toggle(key, *, default: bool, **kwargs):
        def set_value(value):
            set_shared_value(key, value)

        toggle = widgets.ToggleButton(**kwargs, value=get_shared_value(key, default))
        toggle.observe(lambda change: set_value(change["new"]), names="value")
        return toggle

    def interact_state_slider_float(key, *, default: float, **kwargs):
        def set_value(value):
            set_shared_value(key, value)

        return interact(
            set_value,
            value=widgets.FloatSlider(
                **kwargs,
                value=get_shared_value(key, default),
            ),
        )

    def interact_state_slider_int(key, *, default: int, **kwargs):
        def set_value(value):
            set_shared_value(key, value)

        return interact(
            set_value,
            value=widgets.IntSlider(
                **kwargs,
                value=get_shared_value(key, default),
            ),
        )

    def reset_molecule(_):
        app_server.run_command(keys.RESET_COMMAND, {})

    def reset_box(_):
        set_shared_value("scene", [0, 0, 0, 0, 0, 0, 1, 1, 1, 1])

    def copy_box(_):
        nonlocal copied_box
        copied_box = get_shared_value("scene")

    def paste_box(_):
        if copied_box is not None:
            set_shared_value("scene", copied_box)

    def set_force_type(type="spring"):
        set_shared_value("suggested.interaction.type", type)

    def set_simulation_paused(paused=False):
        if not paused:
            app_server.run_command(keys.PLAY_COMMAND, {})
        else:
            app_server.run_command(keys.PAUSE_COMMAND, {})

    def make_switch(i):
        def switch(_):
            app_server.run_command(keys.LOAD_COMMAND, {"index": i})

        return switch

    simulation_names = app_server.run_command(keys.LIST_COMMAND, {}).result()[
        "simulations"
    ]

    switching = []
    for i, name in enumerate(simulation_names):
        button = widgets.Button(description=f"{name}")
        button.on_click(make_switch(i))
        switching.append(button)

    control = [
        make_button(reset_molecule, description="Reset Molecule"),
        make_toggle(set_simulation_paused, description="Pause Simulation"),
        make_state_toggle(
            "suggested.interaction.hydrogens",
            description="Interactable Hydrogens",
            default=True,
        ),
    ]

    box = [
        make_button(copy_box, description="Copy Box"),
        make_button(paste_box, description="Paste Box"),
        make_button(reset_box, description="Reset Box"),
        make_state_toggle(
            "suggested.box.locked", description="Lock Box", default=False
        ),
        make_state_toggle(
            "suggested.box.hidden", description="Hide Box", default=False
        ),
    ]

    interact(
        set_force_type,
        type=widgets.Dropdown(
            options=[
                ("Gaussian", "gaussian"),
                ("Harmonic", "spring"),
                ("Constant", "constant"),
            ],
            description="Force Type",
            value="spring",
        ),
        value=get_shared_value("suggested.interaction.type", "spring"),
    )

    interact_state_slider_int(
        "suggested.interaction.scale",
        default=100,
        description="Force Scale",
        min=1,
        max=1000,
    )
    interact_state_slider_float(
        "suggested.interaction.range",
        default=0.4,
        description="Force Range",
        min=0.1,
        max=2.0,
        step=0.05,
    )
    interact_state_slider_float(
        "suggested.passthrough",
        default=1,
        description="Passthrough",
        min=0,
        max=1,
    )

    display(
        widgets.HBox(
            [
                widgets.VBox(control),
                widgets.VBox(switching),
                widgets.VBox(box),
            ]
        )
    )
