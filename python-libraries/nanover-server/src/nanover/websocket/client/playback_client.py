from nanover.trajectory.frame_server import (
    PLAY_COMMAND_KEY,
    NEXT_COMMAND_KEY,
    STEP_COMMAND_KEY,
    PAUSE_COMMAND_KEY,
    RESET_COMMAND_KEY,
    LOAD_COMMAND_KEY,
    LIST_COMMAND_KEY,
)
from nanover.websocket.client.base_client import WebsocketClient


class PlaybackClient(WebsocketClient):
    """
    Mixin of methods for running playback commands with a WebSocketClient.
    """

    def run_play(self):
        """
        Sends a request to start playing the trajectory to the trajectory service.
        """
        self.run_command_blocking(PLAY_COMMAND_KEY)

    def run_step(self):
        """
        Sends a request to take one step to the trajectory service.
        """
        self.run_command_blocking(STEP_COMMAND_KEY)

    def run_pause(self):
        """
        Sends a request to pause the simulation to the trajectory service.
        """
        self.run_command_blocking(PAUSE_COMMAND_KEY)

    def run_reset(self):
        """
        Sends a request to reset the simulation to the trajectory service.
        """
        self.run_command_blocking(RESET_COMMAND_KEY)

    def run_load(self, index: int):
        """
        Sends a request for the trajectory service to switch to a particular simulation.
        """
        self.run_command_blocking(LOAD_COMMAND_KEY, index=index)

    def run_next(self):
        """
        Sends a request for the trajectory service to switch to the next simulation.
        """
        self.run_command_blocking(NEXT_COMMAND_KEY)

    def run_list(self) -> list[str]:
        """
        Retrieves an ordered list of the available simulations for the load command.
        """
        return self.run_command_blocking(LIST_COMMAND_KEY)["simulations"]
