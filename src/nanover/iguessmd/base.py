from abc import abstractmethod
from os import PathLike

import numpy as np
import numpy.typing as npt


class iGUESSMDSimulation:
    """
    A base class for all iGUESSMD simulations.
    """

    # TODO: Flesh out this base class and move operations common to all
    #  subclasses to this base class

    def __init__(self, name: str | None = None):
        self.name: str | None = name or "Unnamed iGUESSMD simulation"

    @abstractmethod
    def define_iguessmd_simulation_atom_positions_array(self):
        """
        Define the array to which the positions of the atoms with which the
        iGUESSMD force interacts will be saved over the course of the iGUESSMD simulation.
        """

    @abstractmethod
    def add_iguessmd_force_to_system(self):
        """
        Add the required iGUESSMD force to the system, depending on the type of
        iGUESSMD interaction required (single atom or centre-of-mass).
        """

    @abstractmethod
    def run_equilibration_with_initial_restraint(self, n_steps: int):
        """
        Perform an equilibration of the system with the restraint fixed at the
        initial position defined in the iGUESSMD force path.

        :param n_steps: Number of simulation steps to run equilibration for.
        """

    @abstractmethod
    def update_iguessmd_force_position(self):
        """
        Update the position of the iGUESSMD force.
        """

    @abstractmethod
    def run_iguessmd(self, progress_interval: int = 100000):
        """
        Perform an iGUESSMD simulation on the system using the iGUESSMD force and
        path defined, store the positions of the atoms with which the
        iGUESSMD force interacts, and calculate the work done along the reaction
        coordinate.

        :param progress_interval: Interval defining how regularly to print the progress
          of the simulation, in terms of the number of simulation steps.
        """

    @abstractmethod
    def calculate_cumulative_work_done(self):
        """
        Calculate the cumulative work done along the reaction coordinate
        during the iGUESSMD simulation.
        """

    @abstractmethod
    def save_iguessmd_simulation_data(
        self,
        path: PathLike | str | None = None,
        every_nth: int | None = None,
        include_end: bool = False,
        save_work_done: bool = True,
        save_atom_positions: bool = True,
        save_com_positions: bool = True,
        work_done_dtype: npt.DTypeLike = np.float32,
        atom_positions_dtype: npt.DTypeLike = np.float32,
        com_positions_dtype: npt.DTypeLike = np.float32,
    ):
        """
        Save the data produced by the iGUESSMD simulation in binary form that can be read
        into NumPy arrays. The following data are saved, in the order listed below:

        - Trajectories of the atoms to which the iGUESSMD force was applied, in nm
        - Work done along the reaction coordinate defined by the path of the iGUESSMD force, in kJ mol-1

        :param path: Path to the file to which the data will be saved
        :param every_nth:  (int | None) only save the values of work and/or positions at every
          nth point along the trajectory
        :param include_end: (Bool) whether to save final values of the work and/or positions
          regardless of stride defined by every nth
        :param save_work_done: Bool determining whether to save the work done
        :param save_atom_positions: Bool determining whether to save the positions of the atom(s)
        :param save_com_positions: Bool determining whether to save the positions of the COM
        :param work_done_dtype: Data type of the work done array to save
        :param atom_positions_dtype: Data type of the atom positions array to save
        :param com_positions_dtype: Data type of the COM positions array to save
        """

    @abstractmethod
    def save_general_iguessmd_data(self, path: PathLike | str | None = None):
        """
        Saves general data related to the iGUESSMD simulation in binary form that can be read
        into NumPy arrays. The following data are saved, in the order listed below:

        - Indices of the atoms to which the iGUESSMD force is applied
        - Positions defining the path of the iGUESSMD force, in nm
        - Force constant of the iGUESSMD force, in kJ mol-1 nm-2
        - Temperature of the simulation, in Kelvin
        - Time step of the simulation, in picoseconds

        :param path: Path to the file to which the data will be saved.
        """
