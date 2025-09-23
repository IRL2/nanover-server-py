"""
Check a NanoVer recording for errors.

Example when used as a cli:

.. code:: bash

    # show the help
    nanover-check-recording --help

    # check a recording
    nanover-check-recording recording.traj recording.state

"""

import argparse
from pathlib import Path

from nanover.recording2.convert import convert_old_recording


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        dest="paths",
        nargs="+",
        metavar="PATH",
        help="Recording files to convert (one or both of .traj and .state)",
    )
    args = parser.parse_args()

    paths = [Path(path) for path in args.paths]
    traj_path = next((path for path in paths if path.suffix == ".traj"), None)
    state_path = next((path for path in paths if path.suffix == ".state"), None)
    out_path = f"{(traj_path or state_path or "FAILED")}.nanover.zip"

    print(f"Writing recording to {out_path}")
    convert_old_recording(out_path, traj=traj_path, state=state_path)


if __name__ == "__main__":
    main()
