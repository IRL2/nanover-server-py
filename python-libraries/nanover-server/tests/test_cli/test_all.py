import subprocess

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
