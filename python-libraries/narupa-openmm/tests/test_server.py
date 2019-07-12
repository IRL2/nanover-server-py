"""
Tests for the :class:`narupa.openmm.Server` facility.
"""
# Pylint does not recognize pytest fixtures which creates fake warnings.
# pylint: disable=redefined-outer-name,unused-import
# It is expected to access "private" attributes during the tests.
# pylint: disable=protected-access
# Inherited test methods loose the staticmethod decorator. Test method that
# will not be overwritten therefore cannot be staticmethods, even if they do
# not use self.
# pylint: disable=no-self-use
import pytest

from narupa.openmm import Server, NarupaReporter

from simulation_utils import (
    DoNothingReporter,
    basic_simulation,
    serialized_simulation_path,
)
from test_runner import TestRunner


class TestServer(TestRunner):
    """
    Tests for the :class:`Server` class.

    This runs the tests for the :class:`narupa.openmm.Runner` class applied to
    :class:`Server`, as well as the :class:`Server` specific tests.

    When writing a test that requires a port, use the :attr:`next_port` property.
    This property returns the next port starting from :attr:`starting_port` + 1.
    This allows to avoid using twice the same port in two different tests as freeing
    a port takes time. Use the :attr:`last_port` property to get the last port that
    was provided.
    """
    __test__ = True

    _current_port = 6000
    expected_number_of_reporters_verbosity = {
        True: 3,
        False: 2,
    }

    @pytest.fixture
    def runner(self, basic_simulation):
        """
        Setup a :class:`Server` on a basic simulation.

        The simulation has a reporter attached to it to assure removing a reporter
        only removes only that reporter.
        """
        server = Server(basic_simulation, address='[::]', port=self.next_port)
        server.simulation.reporters.append(DoNothingReporter())
        yield server
        server.close()
    
    @property
    def next_port(self):
        # Potentially, a new instance of that class is created for each test.
        # Therefore, we need the port to be stored at the class level, and not
        # at the instance level.
        self.__class__._current_port += 1
        return self._current_port
    
    @property
    def last_port(self):
        return self._current_port

    def test_class(self, runner):
        assert isinstance(runner, Server)

    def test_from_xml_input(self, serialized_simulation_path):
        server = Server.from_xml_input(
            serialized_simulation_path,
            address='[::]', port=self.next_port,
            publish_interval=2,
        )
        n_atoms_in_system = 8
        assert server.simulation.system.getNumParticles() == n_atoms_in_system
        assert isinstance(server._frame_reporter, NarupaReporter)

    def test_default_publishing_frames(self, runner):
        """
        By default, frames are published.
        """
        assert runner.publishing_frames

    def test_make_publish_frames(self, runner):
        """
        :meth:`Server.make_publish_frames` attaches the reporter.
        """
        reporters = runner.simulation.reporters
        runner.publishing_frames = False
        assert not runner.publishing_frames
        runner.make_publish_frames()
        assert runner.publishing_frames
        assert len(reporters) == self.expected_number_of_reporters_verbosity[runner.verbose]

    def test_make_not_publish_frames(self, runner):
        """
        :meth:`Server.make_not_publish_frames` removes the reporter.
        """
        reporters = runner.simulation.reporters
        assert runner.publishing_frames
        runner.make_not_publish_frames()
        assert not runner.publishing_frames
        assert len(reporters) == self.expected_number_of_reporters_verbosity[runner.verbose] - 1

    @pytest.mark.parametrize('initial_value', (True, False))
    @pytest.mark.parametrize('set_value_to', (True, False))
    def test_set_publishing_frames_from_property(self, runner, initial_value, set_value_to):
        """
        The :attr:`Server.publishing_frames` property can attach or detach the
        reporter.
        """
        reporters = runner.simulation.reporters
        base_number_of_reporters = self.expected_number_of_reporters_verbosity[runner.verbose] - 1
        runner.publishing_frames = initial_value
        assert runner.publishing_frames == initial_value
        assert len(reporters) == base_number_of_reporters + int(initial_value)
        runner.publishing_frames = set_value_to
        assert runner.publishing_frames == set_value_to
        assert len(reporters) == base_number_of_reporters + int(set_value_to)

    def test_default_address(self, basic_simulation):
        """
        We can instantiate a server without providing the address.
        """
        # TODO: The address in used should be checked, but I do not know how
        #  to access it. See issue #59.
        server = Server(basic_simulation, port=self.next_port)
        server.close()

    def test_default_port(self, basic_simulation):
        """
        We can instantiate a server without providing the port.
        """
        # TODO: The port in used should be checked, but I do not know how
        #  to access it. See issue #59.
        server = Server(basic_simulation, address='[::]')
        server.close()

    def test_default_host(self, basic_simulation):
        """
        We can instantiate a server without providing neither an address nor a port.
        """
        # TODO: The address and port in used should be checked, but I do not know how
        #  to access them. See issue #59.
        server = Server(basic_simulation)
        server.close()

    def test_context_manager(self, basic_simulation):
        """
        We can use a server as a context manager without an error.
        """
        with Server(basic_simulation, address='[::]', port=self.next_port) as server:
            server.run(n_steps=1)
