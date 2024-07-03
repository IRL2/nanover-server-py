from unittest.mock import mock_open, patch

from nanover.omni.cli import initialise


@patch('builtins.open', mock_open())
def test_record_opens_files():
    with initialise(["--record", "test"]):
        pass

    open.assert_any_call("test.traj", "wb")
    open.assert_any_call("test.state", "wb")
