"""
Tests for :mod:`narupa.openmm.runner`.
"""
# Pylint does not recognize pytest fixtures which creates fake warnings.
# pylint: disable=redefined-outer-name,unused-import
# Inherited test methods loose the staticmethod decorator. Test method that
# will not be overwritten therefore cannot be staticmethods, even if they do
# not use self.
# pylint: disable=no-self-use
import pytest

from narupa.openmm import Runner

from simulation_utils import (
    DoNothingReporter,
    basic_simulation,
    serialized_simulation_path,
)


class TestRunner:
    """
    Tests for the :class:`Runner` class.

    This class can be inherited to test subclasses of :class:`Runner`. The
    :meth:`runner` fixture and the :meth:`test_class` test must then be
    overwritten. If the subclass adds any reporters, then the
    :attr:`expected_expected_number_of_reporters_verbosity`
    class attribute must be overwritten as well to reflect the default number
    of reporters when the verbosity is set to ``True`` or ``False``.
    """
    expected_number_of_reporters_verbosity = {
        True: 2,
        False: 1,
    }

    @pytest.fixture
    def runner(self, basic_simulation):
        """
        Setup a :class:`Runner` on a basic simulation.

        The simulation has a reporter attached to it to assure removing a reporter
        only removes only that reporter.
        """
        runner = Runner(basic_simulation)
        runner.simulation.reporters.append(DoNothingReporter())
        return runner

    @staticmethod
    def test_class(runner):
        """
        Make sure the :meth:`runner` fixture returns an object of the expected
        class.

        This assures that test classes that inherit from that class use their
        own fixture.
        """
        assert isinstance(runner, Runner)

    def test_default_verbosity(self, runner):
        """
        Test that the verbosity is off by default
        """
        assert not runner.verbose

    @pytest.mark.parametrize('initial_value', (True, False))
    @pytest.mark.parametrize('set_value_to', (True, False))
    def test_set_verbosity_from_property(self, runner, initial_value, set_value_to):
        """
        Test that the verbosity can be set from the :attr:`Runner.verbose` property.

        The test makes sure that the value can be set from one value to an other,
        and from one value to itself.
        """
        reporters = runner.simulation.reporters
        runner.verbose = initial_value
        assert runner.verbose == initial_value
        assert len(reporters) == self.expected_number_of_reporters_verbosity[initial_value]
        runner.verbose = set_value_to
        assert runner.verbose == set_value_to
        assert len(reporters) == self.expected_number_of_reporters_verbosity[set_value_to]

    @pytest.mark.parametrize('initial_value', (True, False))
    def test_make_verbose(self, runner, initial_value):
        """
        Test that :meth:`Runner.make_verbose` sets the verbosity on.
        """
        reporters = runner.simulation.reporters
        runner.verbose = initial_value
        assert runner.verbose == initial_value
        runner.make_verbose()
        assert runner.verbose
        assert len(reporters) == self.expected_number_of_reporters_verbosity[True]

    @pytest.mark.parametrize('initial_value', (True, False))
    def test_make_quiet(self, runner, initial_value):
        """
        Test that :meth:`Runner.make_quiet` sets the verbosity off.
        """
        reporters = runner.simulation.reporters
        runner.verbose = initial_value
        assert runner.verbose == initial_value
        runner.make_quiet()
        assert not runner.verbose
        assert len(reporters) == self.expected_number_of_reporters_verbosity[False]

    def test_run(self, runner):
        """
        Test that :meth:`Runner.run` runs the simulation.
        """
        assert runner.simulation.currentStep == 0
        runner.run(n_steps=5)
        assert runner.simulation.currentStep == 5

    def test_from_xml_input(self, serialized_simulation_path):
        """
        Test that a :class:`Runner` can be built from a serialized simulation.
        """
        runner = Runner.from_xml_input(serialized_simulation_path)
        n_atoms_in_system = 8
        assert runner.simulation.system.getNumParticles() == n_atoms_in_system
