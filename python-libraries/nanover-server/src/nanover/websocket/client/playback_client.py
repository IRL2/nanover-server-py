from nanover.trajectory import keys
from nanover.websocket.client.base_client import WebsocketClient


class PlaybackClient(WebsocketClient):
    """
    Mixin of methods for running playback commands with a WebSocketClient.
    """

    def run_play(self):
        """
        Sends a request to start playing the trajectory to the trajectory service.
        """
        self.run_command_blocking(keys.PLAY_COMMAND)

    def run_step(self):
        """
        Sends a request to take one step to the trajectory service.
        """
        self.run_command_blocking(keys.STEP_COMMAND)

    def run_pause(self):
        """
        Sends a request to pause the simulation to the trajectory service.
        """
        self.run_command_blocking(keys.PAUSE_COMMAND)

    def run_reset(self):
        """
        Sends a request to reset the simulation to the trajectory service.
        """
        self.run_command_blocking(keys.RESET_COMMAND)

    def run_load(self, index: int):
        """
        Sends a request for the trajectory service to switch to a particular simulation.
        """
        self.run_command_blocking(keys.LOAD_COMMAND, index=index)

    def run_next(self):
        """
        Sends a request for the trajectory service to switch to the next simulation.
        """
        self.run_command_blocking(keys.NEXT_COMMAND)

    def run_list(self) -> list[str]:
        """
        Retrieves an ordered list of the available simulations for the load command.
        """
        return self.run_command_blocking(keys.LIST_COMMAND)["simulations"]
