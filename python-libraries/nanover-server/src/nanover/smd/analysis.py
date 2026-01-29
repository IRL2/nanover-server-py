import numpy as np

import openmm.unit as unit
from openmm.unit.quantity import Quantity


def boltzmann_constant_in_kJ_mol_K() -> Quantity:
    """
    Returns the Boltzmann constant in units of kJ mol-1 K-1
    """
    kB_J_mol_K = unit.BOLTZMANN_CONSTANT_kB * unit.AVOGADRO_CONSTANT_NA
    return kB_J_mol_K.in_units_of(unit.kilojoule_per_mole / unit.kelvin)


def calculate_beta_mol_kJ(temperature_K: float) -> float:
    """
    Calculate the value of beta (1/(kB*T)) in units of mol kJ-1.

    :param temperature_K: Temperature (in K)
    :return: Beta (1/(kB*T)) in units of mol kJ-1.
    """
    kB_kJ_mol_K = boltzmann_constant_in_kJ_mol_K()
    beta = 1.0 / (kB_kJ_mol_K._value * temperature_K)
    return beta


def _calculate_pmf_second_cumulant(
    work_done_array: np.ndarray, beta: float
) -> np.ndarray:
    """
    Calculate the PMF from the irreversible work done via the second cumulant expansion
    of the Jarzynski equality. This uses the inbuilt Bessel correction in NumPy to unbias
    the calculated variance, assuming that the set of work profiles used to compute the
    PMF is a sample of the (true) larger population.
    WARNING: the irreversible work and beta must have compatible units!
    """
    return np.average(work_done_array, axis=0) - 0.5 * beta * np.var(
        work_done_array, axis=0, ddof=1
    )

def _calculate_pmf_exponential_average(
    work_done_array: np.ndarray, beta: float
):
    """
    Calculate the PMF from the irreversible work done via the exponential average.
    WARNING: the irreversible work and beta must have compatible units!
    """
    return - (1. / beta) * np.log(np.average(np.exp(- beta * work_done_array), axis=0))

def calculate_pmf_second_cumulant_kJ_mol(
    work_done_array_kJ_mol: np.ndarray, temperature_K: float
) -> np.ndarray:
    """
    Calculate the PMF in kJ mol-1 from the irreversible work done (also in kJ mol-1)
    via the second cumulant expansion of the Jarzynski equality.

    :param work_done_array_kJ_mol: N x k array of irreversible work done from a set of
      N SMD simulations performed using a reaction coordinate defined by k positions
      (in kJ mol-1).
    :param temperature_K: Temperature (in K)
    """
    kB_kJ_mol_K = boltzmann_constant_in_kJ_mol_K()
    beta = calculate_beta_mol_kJ(temperature_K)
    return _calculate_pmf_second_cumulant(work_done_array_kJ_mol, beta)

def calculate_pmf_exponential_average_kJ_mol(
    work_done_array_kJ_mol: np.ndarray, temperature_K: float
):
    """
    Calculate the PMF in kJ mol-1 from the irreversible work done (also in kJ mol-1)
    via the exponential average.

    :param work_done_array_kJ_mol: N x k array of irreversible work done from a set of
      N SMD simulations performed using a reaction coordinate defined by k positions
      (in kJ mol-1).
    :param temperature_K: Temperature (in K)
    """
    beta = calculate_beta_mol_kJ(temperature_K)
    return _calculate_pmf_exponential_average(work_done_array_kJ_mol, beta)

def calculate_reaction_coordinate_projections(
    smd_com_coordinates_array: np.ndarray,
    smd_reaction_coordinate: np.ndarray,
    every_nth_point: int | None = None,
) -> np.ndarray:
    """
    Calculate the values of the reaction coordinate for the atom/COM coordinates from
    a set of SMD trajectories, via projection onto the reaction coordinate.
    :param smd_com_coordinates_array: (N x k x 3) array of coordinates defining the trajectories of
      the COM of the group pulled along a k-step SMD reaction coordinate during the N SMD simulations
      (in nm)
    :param smd_reaction_coordinate: (k x 3) array of points defining the SMD reaction coordinate (in nm).
    :param every_nth_point: (int | None) only calculate variance in the value of the reaction coordinate
      at every nth point on the trajectory.
    :return: (N * k) array of projected reaction coordinate values
    """
    # Calculate displacement vectors along SMD reaction coordinate
    displacements = calculate_displacements_along_reaction_coordinate(
        smd_reaction_coordinate
    )

    # Assume "displacement" from final simulated point is equal to the
    # final explicit displacement (restraint has same velocity)
    displacements = np.array([*displacements, displacements[-1]])

    # Calculate normalised displacement vectors
    normalised_displacements = np.array(
        [
            displacements[i] / np.linalg.norm(displacements[i])
            for i in range(displacements.shape[0])
        ]
    )

    # Calculate restraint-atom vectors and reaction coordinate values for each trajectory
    restraint_vectors = (smd_com_coordinates_array - smd_reaction_coordinate)
    if every_nth_point is not None:
        displacements = displacements[::every_nth_point]
        restraint_vectors = restraint_vectors[:, ::every_nth_point]

    i_index_range = restraint_vectors.shape[1]
    if abs(restraint_vectors.shape[1] - normalised_displacements.shape[0]) == 1:
        i_index_range -= 1

    reaction_coordinate_projections = np.array(
        [
            [
                np.dot(restraint_vectors[j, i], normalised_displacements[i])
                for i in range(i_index_range)
            ]
            for j in range(restraint_vectors.shape[0])
        ]
    )

    return reaction_coordinate_projections



def calculate_variance_of_reaction_coordinate(
    smd_com_coordinates_array: np.ndarray,
    smd_reaction_coordinate: np.ndarray,
    every_nth_point: int | None = None,
) -> np.ndarray:
    """
    Calculates the variance in the reaction coordinate value for a set of SMD simulations.

    :param smd_com_coordinates_array: (N x k x 3) array of coordinates defining the trajectories of
      the COM of the group pulled along a k-step SMD reaction coordinate during the N SMD simulations
      (in nm)
    :param smd_reaction_coordinate: (k x 3) array of points defining the SMD reaction coordinate (in nm).
    :param every_nth_point: (int | None) only calculate variance in the value of the reaction coordinate
      at every nth point on the trajectory.
    :return: (N-1) array of the variance in the reaction coordinate as a function of the first N-1
      points of the reaction coordinate (in nm^2)
    """
    reaction_coordinate_projections = (
        calculate_reaction_coordinate_projections(smd_com_coordinates_array, smd_reaction_coordinate, every_nth_point))

    # # Calculate displacement vectors along SMD reaction coordinate
    # displacements = calculate_displacements_along_reaction_coordinate(
    #     smd_reaction_coordinate, every_nth_point=every_nth_point
    # )
    #
    # # Calculate normalised displacement vectors
    # normalised_displacements = np.array(
    #     [
    #         displacements[i] / np.linalg.norm(displacements[i])
    #         for i in range(displacements.shape[0])
    #     ]
    # )
    #
    # # Calculate restraint-atom vectors and reaction coordinate values for each trajectory
    # restraint_vectors = (smd_com_coordinates_array - smd_reaction_coordinate)[:, :-1]
    # if every_nth_point is not None:
    #     restraint_vectors = restraint_vectors[:, ::every_nth_point]
    #
    # i_index_range = restraint_vectors.shape[1]
    # if abs(restraint_vectors.shape[1] - normalised_displacements.shape[0]) == 1:
    #     i_index_range -= 1
    #
    # reaction_coordinate_projections = np.array(
    #     [
    #         [
    #             np.dot(restraint_vectors[j, i], normalised_displacements[i])
    #             for i in range(i_index_range)
    #         ]
    #         for j in range(restraint_vectors.shape[0])
    #     ]
    # )

    return np.var(reaction_coordinate_projections, axis=0)


def calculate_displacements_along_reaction_coordinate(
    smd_reaction_coordinate: np.ndarray,
    every_nth_point: int | None = None,
) -> np.ndarray:
    """
    Calculates the displacements along the reaction coordinate (in nm)

    :param smd_reaction_coordinate: (k x 3) array of points defining the SMD reaction coordinate (in nm).
    :param every_nth_point: (int | None) only return the displacement from one point to the next for
      every nth point of the trajectory
    :return: ((N-1) x 3) array of the displacements between the consecutive points defining the
      SMD reaction coordinate (in nm)
    """
    displacements_along_rc = np.diff(smd_reaction_coordinate, axis=0)
    if every_nth_point is None:
        return displacements_along_rc
    else:
        return displacements_along_rc[::every_nth_point]


def calculate_distance_along_reaction_coordinate(
    smd_reaction_coordinate: np.ndarray,
    every_nth_point: int | None = None,
) -> np.ndarray:
    """
    Calculates the cumulative distance travelled along the reaction coordinate (in nm)

    :param smd_reaction_coordinate: (k x 3) array of points defining the SMD reaction coordinate (in nm).
    :param every_nth_point: (int | None) only return the distance travelled along the reaction coordinate
      for every nth point of the trajectory
    :return: (k) array of distances defining the cumulative distance travelled by the SMD restraint
      along the reaction coordinate (in nm)
    """
    distance_along_rc = np.zeros(smd_reaction_coordinate.shape[0])
    displacements_along_rc = calculate_displacements_along_reaction_coordinate(
        smd_reaction_coordinate
    )
    distance_along_rc[1:] = np.cumsum(np.linalg.norm(displacements_along_rc, axis=1))
    if every_nth_point is None:
        return distance_along_rc
    else:
        return distance_along_rc[::every_nth_point]
