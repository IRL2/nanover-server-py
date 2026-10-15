from abc import abstractmethod

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
        pass

    @abstractmethod
    def check_for_existing_iguessmd_force(self):
        """
        Check whether the loaded simulation already contains an iGUESSMD force of
        the correct type for the simulation type.
        """
        pass

    @abstractmethod
    def add_iguessmd_force_to_system(self):
        """
        Add the required iGUESSMD force to the system, depending on the type of
        iGUESSMD interaction required (single atom or centre-of-mass).
        """
        pass

    @abstractmethod
    def update_iguessmd_force_position(self):
        """
        Update the position of the iGUESSMD force.
        """
        pass

    @abstractmethod
    def calculate_cumulative_work_done(self):
        """
        Calculate the cumulative work done along the reaction coordinate
        during the iGUESSMD simulation.
        """
        pass