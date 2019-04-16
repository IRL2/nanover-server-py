"""
Tests for the :class:`narupa.openmm.Server` facility.
"""
import pytest

from narupa.openmm import Server, NarupaReporter
from narupa.trajectory.frame_server import DEFAULT_PORT

from simulation_utils import (
    DoNothingReporter,
    basic_simulation,
    serialized_simulation_path,
)
from test_runner import TestRunner


class TestServer(TestRunner):
    __test__ = True

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
        server = Server(basic_simulation, address='[::]', port=54321)
        server.simulation.reporters.append(DoNothingReporter())
        return server

    def test_class(self, runner):
        assert isinstance(runner, Server)

    def test_from_xml_input(self, serialized_simulation_path):
        server = Server.from_xml_input(
            serialized_simulation_path,
            address='[::]', port=54321,
            publish_interval=2,
        )
        n_atoms_in_system = 8
        assert server.simulation.system.getNumParticles() == n_atoms_in_system
        assert isinstance(server._frame_reporter, NarupaReporter)

    def test_default_publishing_frames(self, runner):
        assert runner.publishing_frames

    def test_make_publish_frames(self, runner):
        reporters = runner.simulation.reporters
        runner.publishing_frames = False
        assert not runner.publishing_frames
        runner.make_publish_frames()
        assert runner.publishing_frames
        assert len(reporters) == self.expected_number_of_reporters_verbosity[runner.verbose]

    def test_make_not_publish_frames(self, runner):
        reporters = runner.simulation.reporters
        assert runner.publishing_frames
        runner.make_not_publish_frames()
        assert not runner.publishing_frames
        assert len(reporters) == self.expected_number_of_reporters_verbosity[runner.verbose] - 1

    @pytest.mark.parametrize('initial_value', (True, False))
    @pytest.mark.parametrize('set_value_to', (True, False))
    def test_set_publishing_frames_from_property(self, runner, initial_value, set_value_to):
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
        #  to access it.
        try:
            server = Server(basic_simulation, port=8000)
        finally:
            server.close()

    def test_default_port(self, basic_simulation):
        """
        We can instantiate a server without providing the port.
        """
        # TODO: The port in used should be checked, but I do not know how
        #  to access it.
        try:
            server = Server(basic_simulation, address='[::]')
        finally:
            server.close()

    def test_default_host(self, basic_simulation):
        """
        We can instantiate a server without providing neither an address nor a port.
        """
        # TODO: The address and port in used should be checked, but I do not know how
        #  to access them.
        try:
            server = Server(basic_simulation)
        finally:
            server.close()

    def test_context_manager(self, basic_simulation):
        """
        We can use a server as a context manager without an error.
        """
        with Server(basic_simulation, address='[::]', port=8000) as server:
            server.run(n_steps=1)
