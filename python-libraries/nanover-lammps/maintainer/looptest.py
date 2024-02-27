"""
Module containing a simple test program that tests the functionality of the LAMMPS hook,
it allows the maintainer to test the code in an loop in a non embedded python environment
by utilising the mocklammps class.
"""

import time
from nanover.lammps import LammpsImd


def main():
    """
    Test call of the LAMMPS hook routine when running outside of LAMMPS.
    """
    h = LammpsImd()
    print("Starting Trajectory Server")
    while True:
        h.lammps_hook()
        print("FRAME STUFF", h.frame_index, h.frame_data.raw)
        time.sleep(1.0 / 10.0)


if __name__ == "__main__":
    main()
