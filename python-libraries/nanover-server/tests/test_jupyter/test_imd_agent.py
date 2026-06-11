import time
from contextlib import contextmanager
from unittest.mock import patch, ThreadingMock
from dataclasses import dataclass

import mock
import pytest

from nanover.app import OmniRunner
from nanover.imd import ParticleInteraction
from nanover.jupyter import ImdAgent
from nanover.testing import assert_equal_soon
from nanover.testing.asserts import assert_true_soon
from nanover.trajectory import FrameData

FRAMES = [FrameData() for i in range(3)]

for i in range(3):
    FRAMES[i].frame_index = i
    FRAMES[i].particle_count = i

FRAMES[0].kinetic_energy = 9


@dataclass
class Setup:
    imd_runner: OmniRunner
    imd_agent: ImdAgent


@contextmanager
def make_setup():
    with (
        OmniRunner.with_basic_server(port=0) as imd_runner,
        ImdAgent.from_runner(imd_runner) as imd_agent,
    ):
        yield Setup(imd_runner, imd_agent)


@pytest.fixture
def setup():
    with make_setup() as setup:
        yield setup


def test_frames_received(setup):
    """Test that each frame published causes a call of `update_interactions`"""
    with patch.object(
        setup.imd_agent, "update_interactions", new=ThreadingMock()
    ) as update_interactions:
        setup.imd_agent.start()
        for frame in FRAMES:
            setup.imd_runner.app_server.frame_publisher.send_frame(frame)
            full_frame = setup.imd_runner.app_server.frame_publisher.current_frame

            update_interactions.wait_until_called()

            update_interactions.assert_called_with(
                full_frame=full_frame, frame_update=mock.ANY
            )
            update_interactions.reset_mock()


def test_frames_pause(setup):
    """Test frames published when paused do not cause a call of `update_interactions`"""
    with patch.object(
        setup.imd_agent, "update_interactions", new=ThreadingMock()
    ) as update_interactions:
        setup.imd_agent.paused = True
        setup.imd_agent.start()

        for frame in FRAMES:
            setup.imd_runner.app_server.frame_publisher.send_frame(frame)
            time.sleep(0)

        setup.imd_agent.paused = False
        setup.imd_runner.app_server.frame_publisher.send_frame(FRAMES[-1])

        update_interactions.wait_until_called()
        full_frame = setup.imd_runner.app_server.frame_publisher.current_frame
        update_interactions.assert_called_once_with(
            full_frame=full_frame, frame_update=mock.ANY
        )


def test_interaction_cleanup(setup):
    """Test no interactions are left after closing agent"""

    assert len(setup.imd_runner.app_server.imd.active_interactions) == 0

    with setup.imd_agent.interactions as interactions:
        for i in range(32):
            interactions.update_interaction(f"test.{i}", ParticleInteraction())

    assert_equal_soon(
        lambda: len(setup.imd_runner.app_server.imd.active_interactions), lambda: 32
    )

    with setup.imd_agent.interactions as interactions:
        for i in range(16):
            interactions.remove_interaction(f"test.{i}")

    assert_equal_soon(
        lambda: len(setup.imd_runner.app_server.imd.active_interactions), lambda: 16
    )

    setup.imd_agent.close()

    assert_equal_soon(
        lambda: len(setup.imd_runner.app_server.imd.active_interactions), lambda: 0
    )
