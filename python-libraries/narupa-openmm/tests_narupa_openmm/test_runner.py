# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.
"""
Tests for :mod:`narupa.openmm.runner`.
"""
# Pylint does not recognize pytest fixtures which creates fake warnings.
# pylint: disable=redefined-outer-name,unused-import
# Inherited test methods loose the staticmethod decorator. Test method that
# will not be overwritten therefore cannot be staticmethods, even if they do
# not use self.
# pylint: disable=no-self-use
import time
import statistics
import math

import pytest

from narupa.app import NarupaImdClient
from narupa.openmm import (
    OpenMMRunner,
    SET_FRAME_INTERVAL_COMMAND_KEY, GET_FRAME_INTERVAL_COMMAND_KEY,
    SET_FORCE_INTERVAL_COMMAND_KEY, GET_FORCE_INTERVAL_COMMAND_KEY,
)
from narupa.openmm.imd import add_imd_force_to_system
from narupa.trajectory.frame_server import (
    PLAY_COMMAND_KEY,
    PAUSE_COMMAND_KEY,
    RESET_COMMAND_KEY,
    STEP_COMMAND_KEY,
    SET_DYNAMICS_INTERVAL_COMMAND_KEY,
    GET_DYNAMICS_INTERVAL_COMMAND_KEY,
)
from narupa.essd import DiscoveryClient

from .simulation_utils import (
    DoNothingReporter,
    basic_simulation_with_imd_force,
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
        True: 3,
        False: 2,
    }

    def make_runner(self, simulation, name=None):
        runner = OpenMMRunner(simulation, port=0, name=name)
        runner.simulation.reporters.append(DoNothingReporter())
        return runner

    @pytest.fixture
    def runner(self, basic_simulation_with_imd_force):
        """
        Setup a :class:`Runner` on a basic simulation.

        The simulation has a reporter attached to it to assure removing a
        reporter only removes only that reporter.
        """
        simulation, _ = basic_simulation_with_imd_force
        runner = self.make_runner(simulation)
        yield runner
        runner.close()

    @pytest.fixture
    def client_runner(self, runner):
        runner_port = runner.app.port
        with NarupaImdClient.connect_to_single_server(port=runner_port) as client:
            yield client, runner

    @staticmethod
    def test_class(runner):
        """
        Make sure the :meth:`runner` fixture returns an object of the expected
        class.

        This assures that test classes that inherit from that class use their
        own fixture.
        """
        assert isinstance(runner, OpenMMRunner)

    def test_app_deprecated(self, runner):
        assert runner.app is runner.app_server
        with pytest.deprecated_call():
            _ = runner.app

    def test_simulation_without_imd_force(self, basic_simulation):
        """
        When created on a simulation without imd force, the runner fails.
        """
        with pytest.raises(ValueError):
            OpenMMRunner(basic_simulation, port=0)

    def test_simulation_multiple_imd_force(self, caplog, basic_simulation):
        """
        When created on a simulation with more than one imd force, the runner
        issues a warning.
        """
        # The forces added to the system will not be accounted for when running
        # the dynamics until the context is reset as the system is already
        # compiled in a context. It does not matter here, as the force is still
        # listed in the system, which is what we check.
        system = basic_simulation.system
        for _ in range(3):
            add_imd_force_to_system(system)

        runner = OpenMMRunner(basic_simulation, port=0)
        runner.close()
        assert 'More than one force' in caplog.text

    @pytest.mark.serial
    @pytest.mark.parametrize('server_name', ("Server 1", "Server 2"))
    def test_discovery_with_client(self, server_name, basic_simulation_with_imd_force):
        simulation, _ = basic_simulation_with_imd_force
        with self.make_runner(simulation, name=server_name) as runner:
            with DiscoveryClient() as client:
                # There may be servers running already, we only want to look at the
                # one we created in that test. We select it by name.
                servers = set(client.search_for_services(search_time=0.8, interval=0.01))
                relevant_servers = [server for server in servers if server.name == server_name]
                assert len(relevant_servers) == 1

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
        expected_number = self.expected_number_of_reporters_verbosity[initial_value]
        assert len(reporters) == expected_number
        runner.verbose = set_value_to
        assert runner.verbose == set_value_to
        expected_number = self.expected_number_of_reporters_verbosity[set_value_to]
        assert len(reporters) == expected_number

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
        runner.run(steps=5)
        assert runner.simulation.currentStep == 5

    def test_from_xml_input(self, serialized_simulation_path):
        """
        Test that a :class:`Runner` can be built from a serialized simulation.
        """
        n_atoms_in_system = 8
        with OpenMMRunner.from_xml_input(serialized_simulation_path, port=0) as runner:
            assert runner.simulation.system.getNumParticles() == n_atoms_in_system

    @pytest.mark.parametrize('name, target_attribute', (
        ('frame_interval', 'frame_interval'),
        ('force_interval', 'force_interval'),
    ))
    def test_interval_get(self, runner, name,  target_attribute):
        """
        The shortcut the the NarupaImdReporter intervals return the expected
        values.
        """
        attribute = getattr(runner.reporter, target_attribute)
        assert getattr(runner, name) == attribute

    @pytest.mark.parametrize('name, target_attribute', (
            ('frame_interval', 'frame_interval'),
            ('force_interval', 'force_interval'),
    ))
    def test_interval_set(self, runner, name, target_attribute):
        """
        The shortcut the the NarupaImdReporter intervals set the expected
        values.
        """
        setattr(runner.reporter, name, 70)
        assert getattr(runner.reporter, name) == 70

    @pytest.mark.parametrize('is_verbose, expectation', (
            (True, 10),
            (False, 0),
    ))
    def test_verbosity_interval_get(self, runner, is_verbose, expectation):
        """
        The shortcut to get the interval for the verbosity print is correct.
        """
        runner.verbose = is_verbose
        assert runner.verbosity_interval == expectation

    @pytest.mark.parametrize('is_verbose', (True, False))
    @pytest.mark.parametrize('interval', (3, 70))
    def test_verbosity_interval_set_non_zero(self, runner, interval, is_verbose):
        """
        Setting the verbosity interval to a non-zero value with the runner
        shortcut alters the reporter and sets the verbosity on.
        """
        runner.verbose = is_verbose
        runner.verbosity_interval = interval
        assert runner._verbose_reporter._reportInterval == interval
        assert runner.verbose is True

    @pytest.mark.parametrize('is_verbose', (True, False))
    def test_verbosity_interval_set_zero(self, runner, is_verbose):
        """
        Setting the verbosity interval to 0 with the runner shortcut alters
        the reporter and sets the verbosity off.
        """
        runner.verbose = is_verbose
        runner.verbosity_interval = 0
        assert runner.verbose is False

    def test_run_non_blocking(self, runner):
        """
        The runner can be run in the background.
        """
        runner.run(100, block=False)
        # Here we count on the context switching to the assertions before
        # finishing the run.
        assert runner.simulation.currentStep < 100
        assert runner._run_task is not None
        assert not runner._cancelled

    def test_cancel_run(self, runner):
        """
        A runner running in the background can be stopped.
        """
        runner.run(block=False)
        assert runner.is_running
        runner.cancel_run(wait=True)
        assert runner.is_running is False

    def test_cancel_never_running(self, runner):
        """
        Cancelling the run of a non-running runner does not raise an error.
        """
        # Cancelling a non running simulation should not raise an exception.
        runner.cancel_run()

    def test_run_twice(self, runner):
        """
        Calling the run method on a running runner raises an exception.
        """
        runner.run(block=False)
        with pytest.raises(RuntimeError):
            runner.run(block=False)

    def test_step(self, runner):
        """
        Calling the step method advances the simulation by the frame_interval.
        """
        step_count = runner.simulation.currentStep
        runner.step()
        assert runner.simulation.currentStep == step_count + runner.frame_interval

    def test_multiple_steps(self, runner):
        """
        The step method can be called multiple times.
        """
        num_steps = 10
        runner.run(block=False)
        runner.cancel_run(wait=True)
        step_count = runner.simulation.currentStep
        for i in range(num_steps):
            runner.step()
        expected_step = step_count + (num_steps * runner.frame_interval)
        assert runner.simulation.currentStep == expected_step

    def test_pause(self, runner):
        """
        The pause method actually pauses the simulation.
        """
        runner.run(block=False)
        runner.pause()
        step_count = runner.simulation.currentStep
        assert runner._run_task.done() is True
        time.sleep(0.1)
        assert runner.simulation.currentStep == step_count

    def test_play(self, runner):
        """
        The play method restart the simulation after a pause.
        """
        runner.run(block=False)
        assert runner.is_running
        runner.pause()
        assert runner.is_running is False
        runner.play()
        assert runner.is_running

    def test_play_twice(self, runner):
        """
        The play method can be called multiple times in succession.
        """
        runner.play()
        assert runner.is_running
        runner.play()
        assert runner.is_running

    @pytest.mark.timeout(1)
    def test_play_command(self, client_runner):
        """
        The play command starts the simulation.
        """
        client, runner = client_runner
        assert not runner.is_running
        client.run_command(PLAY_COMMAND_KEY)
        while True:
            if runner.is_running:
                return

    @pytest.mark.timeout(1)
    def test_pause_command(self, client_runner):
        """
        The pause commands pauses the simulation.
        """
        client, runner = client_runner
        runner.run()
        client.run_command(PAUSE_COMMAND_KEY)
        while runner.is_running:
            continue
        step_count = runner.simulation.currentStep
        time.sleep(0.1)
        assert runner.simulation.currentStep == step_count

    @pytest.mark.timeout(1)
    @pytest.mark.parametrize('running_before', (True, False))
    def test_reset_command(self, client_runner, running_before):
        """
        The reset command calls the reset method and restores the playing status.
        """
        client, runner = client_runner
        runner.run(10)
        if running_before:
            runner.run()

        reset = False

        def on_reset():
            nonlocal reset
            reset = True

        runner.on_reset.add_callback(on_reset)
        client.run_command(RESET_COMMAND_KEY)

        while not reset:
            continue
        assert runner.is_running == running_before

    @pytest.mark.timeout(1)
    def test_step_command(self, client_runner):
        """
        The step command advances and pauses the simulation.
        """
        client, runner = client_runner
        step_count = runner.simulation.currentStep
        client.run_command(STEP_COMMAND_KEY)
        time.sleep(0.1)
        assert runner.is_running is False
        assert runner.simulation.currentStep == step_count + runner.frame_interval

    @pytest.mark.parametrize('target_interval', (5, 20))
    @pytest.mark.parametrize('setter, command', (
            ('frame_interval', GET_FRAME_INTERVAL_COMMAND_KEY),
            ('force_interval', GET_FORCE_INTERVAL_COMMAND_KEY),
    ))
    def test_get_interval_command(
            self, client_runner, target_interval, setter, command):
        client, runner = client_runner
        setattr(runner, setter, target_interval)
        result = client.run_command(command)
        assert result == {'interval': target_interval}

    @pytest.mark.parametrize('target_interval', (5, 20))
    @pytest.mark.parametrize('getter, command', (
            ('frame_interval', SET_FRAME_INTERVAL_COMMAND_KEY),
            ('force_interval', SET_FORCE_INTERVAL_COMMAND_KEY),
    ))
    def test_set_interval_command(
            self, client_runner, target_interval, getter, command):
        client, runner = client_runner
        client.run_command(command, interval=target_interval)
        assert getattr(runner, getter) == target_interval

    def test_get_dynamics_interval(self, runner):
        assert runner.dynamics_interval == runner._variable_interval_generator.interval

    def test_set_dynamics_interval(self, runner):
        value = 20.1
        runner.dynamics_interval = value
        assert runner._variable_interval_generator.interval == pytest.approx(value)

    def test_get_dynamics_interval_command(self, client_runner):
        client, runner = client_runner
        result = client.run_command(GET_DYNAMICS_INTERVAL_COMMAND_KEY)
        time.sleep(0.1)
        assert result == {"interval": pytest.approx(runner.dynamics_interval)}

    def test_set_dynamics_interval_command(self, client_runner):
        value = 10.2
        client, runner = client_runner
        client.run_command(SET_DYNAMICS_INTERVAL_COMMAND_KEY, interval=value)
        time.sleep(0.1)
        assert runner.dynamics_interval == pytest.approx(value)

    @pytest.mark.parametrize("fps", (1, 5, 10, 30))
    @pytest.mark.parametrize("frame_interval", (1, 5, 10))
    def test_throttling(self, client_runner, fps, frame_interval):
        """
        The runner uses the requested MD throttling.

        Here we make sure the runner throttles the dynamics according to the
        dynamics interval. However, we only guarantee that the target dynamics
        interval is a minimum (the MD engine may not be able to produce frames
        fast enough), also we accept some leeway.
        """
        duration = 0.5
        dynamics_interval = 1 / fps
        client, runner = client_runner
        runner.dynamics_interval = dynamics_interval
        runner.frame_interval = frame_interval

        # The frame interval is only taken into account at the end of the
        # current batch of frames. Here we produce one batch of frame and
        # only after that subscribe to the frames.
        runner.run(steps=frame_interval, block=True)
        client.subscribe_to_all_frames()

        runner.run()
        time.sleep(duration)
        runner.cancel_run(wait=True)

        timestamps = [frame.server_timestamp for frame in client.frames]
        deltas = [
            timestamps[i] - timestamps[i - 1]
            for i in range(1, len(timestamps))
        ]
        # The interval is not very accurate. We only check that the observed
        # interval is greater than the expected one and we accept some
        # deviation.
        assert all(delta >= dynamics_interval * 0.90 for delta in deltas)
