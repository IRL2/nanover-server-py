import subprocess
from pathlib import Path

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


@pytest.mark.skip(reason="seems to hang")
def test_split_recording(tmp_path_factory):
    prev_path = (
        Path(__file__).parent.parent
        / "test_omni"
        / "nanotube-example-recording.nanover.zip"
    )
    next_path = tmp_path_factory.mktemp("recording") / "recording.nanover.zip"

    with open(prev_path, "rb") as prev_file, open(next_path, "wb") as next_file:
        next_file.write(prev_file.read())

    assert (
        subprocess.run(
            ["nanover-split-recording", next_path, "-n", "frame.index"]
        ).returncode
        == 0
    )
