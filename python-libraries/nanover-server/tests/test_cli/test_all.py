import subprocess
from pathlib import Path

from nanover.app.cli.server_cli import initialise_runner, handle_user_arguments

import pytest

CLIS = [
    "nanover-server",
    "nanover-record",
    "nanover-essd-list",
    "nanover-dump-state",
    "nanover-split-recording",
]


@pytest.mark.parametrize("cli", CLIS)
def test_cli_help(cli):
    assert subprocess.run([cli, "--help"]).returncode == 0


def test_split_recording(tmp_path_factory):
    prev_path = (
        Path(__file__).parent.parent
        / "test_omni"
        / "nanotube-example-recording.nanover.zip"
    )
    next_path = tmp_path_factory.mktemp("recording") / "recording.nanover.zip"

    with (
        open(prev_path, "rb") as prev_file,
        open(next_path, "wb") as next_file,
    ):
        next_file.write(prev_file.read())

    assert (
        subprocess.run(
            ["nanover-split-recording", next_path, "-n", "frame.index"]
        ).returncode
        == 0
    )


def test_load_mda_structure():
    """Check that static structures can be loaded from commandline."""
    pdb_path = Path(__file__).parent.parent / "test_mdanalysis/17-ala.pdb"
    assert pdb_path.exists(), "Path is incorrect, this file should exist."

    args = handle_user_arguments(args=["--mda", str(pdb_path)])

    print(args)
    with initialise_runner(args) as runner:
        assert len(runner.simulations) == 1
        runner.close()


def test_load_mda_trajectory():
    """Check that trajectories can be loaded from commandline."""
    topology_path = (
        Path(__file__).parent.parent.parent.parent.parent
        / "tutorials/mdanalysis/files/3TI6_ose_wt.pdb"
    )
    trajectory_path = topology_path.with_name("ose_wt.dcd")
    assert topology_path.exists(), "Path is incorrect, this file should exist."
    assert trajectory_path.exists(), "Path is incorrect, this file should exist."

    args = handle_user_arguments(
        args=["--mda", str(topology_path), "--mda-traj", str(trajectory_path)]
    )

    with initialise_runner(args) as runner:
        assert len(runner.simulations) == 1
        assert runner.simulations[0].universe.trajectory.n_frames == 24
        runner.close()


def test_load_mulitple_coordinate_files():
    topology_path = (
        Path(__file__).parent.parent.parent.parent.parent
        / "tutorials/mdanalysis/files/3TI6_ose_wt.pdb"
    )
    trajectory_path = topology_path.with_name("ose_wt.dcd")
    assert topology_path.exists(), "Path is incorrect, this file should exist."
    assert trajectory_path.exists(), "Path is incorrect, this file should exist."

    args = handle_user_arguments(
        args=[
            "--mda",
            str(topology_path),
            "--mda-traj",
            str(trajectory_path),
            str(trajectory_path),
        ]
    )

    with initialise_runner(args) as runner:
        assert len(runner.simulations) == 1
        assert runner.simulations[0].universe.trajectory.n_frames == 48
        runner.close()


def test_load_multiple_mda_trajectories():
    """Checks that mixed sets of trajectories and static structures can be loaded."""
    pdb_path = Path(__file__).parent.parent / "test_mdanalysis/17-ala.pdb"
    topology_path = (
        Path(__file__).parent.parent.parent.parent.parent
        / "tutorials/mdanalysis/files/3TI6_ose_wt.pdb"
    )
    trajectory_path = topology_path.with_name("ose_wt.dcd")

    assert topology_path.exists(), "Path is incorrect, this file should exist."
    assert trajectory_path.exists(), "Path is incorrect, this file should exist."
    assert pdb_path.exists(), "Path is incorrect, this file should exist."

    args = handle_user_arguments(
        args=[
            "--mda",
            str(topology_path),
            "--mda-traj",
            str(trajectory_path),
            "--mda",
            str(pdb_path),
        ]
    )

    with initialise_runner(args) as runner:
        assert len(runner.simulations) == 2
        assert runner.simulations[0].name == str(topology_path)
        assert runner.simulations[1].name == str(pdb_path)
        runner.close()
