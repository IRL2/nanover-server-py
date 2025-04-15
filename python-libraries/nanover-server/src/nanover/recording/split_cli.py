import argparse
from pathlib import Path
from typing import Iterable

from nanover.recording.reading import iter_full_view_max
from nanover.trajectory.frame_data import SIMULATION_COUNTER


def handle_user_arguments(args=None):
    """
    Parse the arguments from the command line.

    :return: The namespace of arguments read from the command line.
    """
    parser = argparse.ArgumentParser(description="Split long recordings into multiple recordings.")

    parser.add_argument(
        "paths",
        nargs="+",
    )

    arguments = parser.parse_args(args)
    return arguments


def iter_prev_next(seq):
    prev_item = next(seq)
    for next_item in seq:
        yield prev_item, next_item
        prev_item = next_item


if __name__ == "__main__":
    arguments = handle_user_arguments()

    paths = [Path(path) for path in arguments.paths]
    traj_path = next((path for path in paths if path.suffix == ".traj"), None)
    state_path = next((path for path in paths if path.suffix == ".state"), None)

    sequences = []
    prev_updates = []

    for prev_data, next_data in iter_prev_next(iter_full_view_max(traj=traj_path, state=state_path)):
        prev_updates.append(prev_data)

        if next_data.frame and SIMULATION_COUNTER in next_data.frame.values:
            sequences.append(prev_updates)
            prev_updates = []
            print("DIVIDE", next_data.frame.values[SIMULATION_COUNTER])
