import pytest
from MDAnalysis import Universe

from nanover.mdanalysis import UniverseSimulation

from .common import make_loaded_sim, PDB_PATH, DCD_PATH


@pytest.fixture
def example_playback_mda():
    with make_loaded_sim(make_example_playback_mda()) as sim:
        yield sim


def make_example_playback_mda():
    universe = Universe(PDB_PATH, DCD_PATH)
    return UniverseSimulation.from_universe(universe)


def test_step_all_frames(example_playback_mda):
    """
    Test that stepping until frame loops sees all frames.
    """
    trajectory = example_playback_mda.universe.trajectory
    prev_frame, next_frame = -1, 0

    frames = []

    while next_frame >= prev_frame:
        prev_frame = trajectory.ts.frame
        example_playback_mda.advance_by_one_step()
        next_frame = trajectory.ts.frame
        frames.append(prev_frame)

    assert frames == list(range(len(trajectory)))


def test_loop_all_frames(example_playback_mda):
    """
    Test that advancing by less than dt until frame loops sees all frames.
    """
    trajectory = example_playback_mda.universe.trajectory
    prev_frame, next_frame = -1, 0
    dt = trajectory.dt * 0.9

    frames = []

    while next_frame >= prev_frame:
        prev_frame = trajectory.ts.frame
        example_playback_mda.advance_by_seconds(dt)
        next_frame = trajectory.ts.frame

        if next_frame != prev_frame:
            frames.append(prev_frame)

    assert frames == list(range(len(trajectory)))


def test_last_frame_duration(example_playback_mda):
    """
    Test that last frame lasts roughly the right amount of time.
    """
    trajectory = example_playback_mda.universe.trajectory

    # skip to last frame
    for _ in range(len(trajectory) - 1):
        example_playback_mda.advance_by_one_step()
    assert trajectory.ts.frame == len(trajectory) - 1

    example_playback_mda.advance_by_seconds(trajectory.dt * 0.99)
    assert trajectory.ts.frame == len(trajectory) - 1
