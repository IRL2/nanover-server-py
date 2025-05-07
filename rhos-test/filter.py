"""
Reads the trajectory time series as recorded by nanover, removes all instances of the simulation counter key, and writes
 the output to a new file.

Example usage:

.. code:: bash

    # show the help
    python remove_sim_counter_from_trajectory.py --help

    # remove the simulation counter key and output to a new file
    python remove_sim_counter_from_trajectory.py my_file.traj

"""

from nanover.recording.reading import iter_recording_entries
from nanover.recording.writing import record_entries
from nanover.trajectory.frame_data import SIMULATION_COUNTER
from nanover.protocol.trajectory import GetFrameResponse
import argparse


def remove_simulation_counter_from_frame(frame_response):
    """
    Remove the simulation counter key from the frame if it is present.
    """
    # If the simulation counter key exists, delete it
    if SIMULATION_COUNTER in frame_response.frame.values:
        del frame_response.frame.values[SIMULATION_COUNTER]

    return frame_response


def strip_simulation_counter_from_recording(recording_path_traj):
    """
    Remove any instances of the simulation counter key from the recorded trajectory and output to a new file with ending
     '-stripped'.
    """
    INPATH = recording_path_traj
    OUTPATH = INPATH.replace(".traj", "-stripped.traj")

    with open(INPATH, "rb") as infile, open(OUTPATH, "wb") as outfile:
        entries = iter_recording_entries(infile, GetFrameResponse)
        stripped = ((timestamp, remove_simulation_counter_from_frame(frame_response)) for timestamp, frame_response in entries)
        record_entries(outfile, stripped)


def command_line():
    """
    Run this script from the command line, giving the path to the trajectory recording as an argument.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('path', help='Path to the trajectory recording file to read.')
    args = parser.parse_args()
    strip_simulation_counter_from_recording(args.path)


if __name__ == '__main__':
    command_line()