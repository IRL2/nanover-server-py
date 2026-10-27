from nanover.imd import ParticleInteraction


class Mode:
    def on_enter(self):
        pass

    def on_exit(self):
        pass

    def on_button_pressed(self, *, key: str, cursor: dict, button: str):
        pass

    def on_button_released(self, *, key: str, cursor: dict, button: str):
        pass

    def on_cursor_updated(self, *, key: str, cursor: dict):
        pass

    def on_interaction_started(self, *, key: str, interaction: ParticleInteraction):
        pass

    def on_interaction_stopped(self, *, key: str, interaction: ParticleInteraction):
        pass

    def on_interaction_updated(self, *, key: str, interaction: ParticleInteraction):
        pass
