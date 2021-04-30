from abc import ABCMeta, abstractmethod
from typing import Optional


class NarupaRunner(metaclass=ABCMeta):

    @property
    def address(self):
        """
        Gets the URL or IP address the server is running at.
        """
        return self.app_server.address

    @property
    def port(self):
        """
        Gets the port the server is running on.
        """
        return self.app_server.port

    @property
    def name(self):
        """
        Get the name of the server.
        """
        return self.app_server.name

    @property
    @abstractmethod
    def app_server(self):
        """Get the underlying app of the runner."""
        raise NotImplementedError()

    @abstractmethod
    def run(
            self,
            steps: Optional[int] = None,
            block: Optional[bool] = None,
    ) -> None:
        """
        Runs the molecular dynamics.

        :param steps: If passed, will run the given number of steps, otherwise
            will run forever on a background thread and immediately return.
        :param block: If ``False``, run in a separate thread. By default, "block"
            is ``None``, which means it is automatically set to ``True`` if a
            number of steps is provided and to ``False`` otherwise.
        """
        raise NotImplementedError()

    @abstractmethod
    def step(self):
        """
        Take a single step of the simulation and stop.

        This method is called whenever a client runs the step command,
        described in :mod:narupa.trajectory.frame_server.
        """
        raise NotImplementedError()

    @abstractmethod
    def pause(self):
        """
        Pause the simulation, by cancelling any current run.

        This method is called whenever a client runs the pause command,
        described in :mod:narupa.trajectory.frame_server.
        """
        raise NotImplementedError()

    @abstractmethod
    def play(self):
        """
        Run the simulation indefinitely

        Cancels any current run and then begins running the simulation on a background thread.

        This method is called whenever a client runs the play command,
        described in :mod:narupa.trajectory.frame_server.
        """
        raise NotImplementedError()

    @abstractmethod
    def reset(self):
        raise NotImplementedError()

    @abstractmethod
    def close(self):
        """
        Closes the connection and stops the dynamics.
        """
        raise NotImplementedError()

    @property
    @abstractmethod
    def is_running(self):
        """
        Whether or not the molecular dynamics is currently running on a background thread or not.

        :return: `True`, if molecular dynamics is running, `False` otherwise.
        """
        raise NotImplementedError()

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        self.close()
