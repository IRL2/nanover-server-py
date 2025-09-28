"""
Read the shared state time series as recorded by NanoVer.

Example when used as a cli:

.. code:: bash

    # show the help
    nanover-dump-state --help

    # show the timestamps and state updates one record per line
    nanover-dump-state my_file.state
    # show the timestamps and the state updates in human-readable form
    nanover-dump-state --pretty my_file.state

    # show the timestamps and full states one record per line
    nanover-dump-state --full my_file.state
    # show the timestamps and full state in human-readable form
    nanover-dump-state --full --pretty my_file.state

"""

import argparse
from pprint import pprint

from nanover.recording2.reading import NanoverRecordingReader


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--full",
        action="store_true",
        default=False,
        help="Display the aggregated state instead of the state updates.",
    )
    parser.add_argument(
        "--pretty",
        action="store_true",
        default=False,
        help="Display the state in a more human-readable way.",
    )
    parser.add_argument("path", help="Path to the file to read.")
    args = parser.parse_args()

    reader = NanoverRecordingReader.from_path(args.path)

    if args.full:
        stream = (
            (entry.metadata["timestamp"], state)
            for entry, state in reader.iter_states_full()
        )
    else:
        stream = (
            (entry.metadata["timestamp"], update.updates)
            for entry, update in reader.iter_state_updates()
        )

    for timestamp, state in stream:
        if args.pretty:
            print(f"---- {timestamp} ---------")
            pprint(state)
        else:
            print(timestamp, state)


if __name__ == "__main__":
    main()
