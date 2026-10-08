"""
Tests for :mod:`nanover.iguessmd.analysis`.

Things to test for analysis functions:
- General iGUESSMD data loads correctly [ ]
- Boltzmann constant is calculated correctly in units of kJ mol-1 K-1 [√]
- Beta (1 / (kB * T)) is correctly calculated in units of mol kJ-1 [√]
- Calculation of PMF via exponential average returns expected result [√]
- Calculation of PMF via second cumulant approximation returns expected result [√]
- Calculation of reaction coordinate projections works as expected [√]
- Calculation of displacements along reaction path works as expected [√]
- Calculation of distance along reaction coordinate works as expected [√]
"""

import pytest
import numpy as np

from itertools import product

from nanover.iguessmd.analysis import *

from .iguessmd_test_utilities import define_circular_path
from scipy.special import logsumexp

KB_KJ_MOL_K_VALUE = 0.008314462618

TEST_iGUESSMD_LINEAR_PATH = np.array(
    [np.linspace(0.05, 1.05, 101), np.zeros(101), np.zeros(101)]
).transpose()
TEST_iGUESSMD_CIRCULAR_PATH, TEST_iGUESSMD_CIRCULAR_PATH_TANGENTS = define_circular_path(10)

TEST_iGUESSMD_POSITION_SHIFTS = [
    np.array([0.0, 0.0, 0.0]),
    np.array([1.0, 0.0, 0.0]),
    np.array([0.0, 1.0, 0.0]),
    np.array([0.0, 0.0, 1.0]),
    np.array([2.0, 0.0, 0.0]),
    np.array([-1.75, -3.0, 5.263]),
]


def test_calculate_boltzmann_constant_in_kJ_mol_K():
    """
    Check that the Boltzmann constant is correctly calculated in units of kJ mol-1 K-1
    """
    kB_calculated = boltzmann_constant_in_kJ_mol_K()
    assert kB_calculated._value == pytest.approx(KB_KJ_MOL_K_VALUE, rel=1e-10)


@pytest.mark.parametrize(
    "temperature_K", [-1.0, 0.0, 1, np.pi, 273.15, 300.0, 12345.6789]
)
def test_calculate_beta_mol_kJ(temperature_K):
    """
    Check that the value of (kB * T)^{-1} is correctly calculated in units
    of mol kJ-1 for positive (non-negative and non-zero) temperatures.
    """
    # Check that an error is thrown for negative or zero temperatures
    if temperature_K <= 0.0:
        with pytest.raises(ValueError):
            calculate_beta_mol_kJ(temperature_K)
        return

    beta_calculated = calculate_beta_mol_kJ(temperature_K)
    beta_expected = 1.0 / (KB_KJ_MOL_K_VALUE * temperature_K)
    assert beta_calculated == pytest.approx(beta_expected, rel=1e-10)


@pytest.mark.parametrize(
    "position_shifts, n_traj",
    product(TEST_iGUESSMD_POSITION_SHIFTS, [1, 15]),
)
def test_calculate_reaction_coordinate_projections(position_shifts, n_traj):
    """
    Check that the reaction coordinate projections are calculated correctly, both in the
    case that the user returns the projections of the full trajectory and in the case
    that the user specifies a non-zero every_nth. Tests written for a linear reaction path and
    a set of positions shifted from the reaction path by a constant value
    (position_shifts) from n_traj trajectories arranged in an (n_traj x k x 3) array.
    """
    # Define every_nth and length of positions array for each trajectory
    every_nth = 11
    array_len = TEST_iGUESSMD_LINEAR_PATH.shape[0]

    # Define positions array of n_traj trajectories
    iguessmd_positions = np.array(
        [TEST_iGUESSMD_LINEAR_PATH + position_shifts for i in range(n_traj)]
    )

    # Calculate expected projections array
    expected_projections = np.zeros((n_traj, array_len))
    expected_projections[:, :] = position_shifts[0]

    # Calculate full projections array
    calculated_projections = calculate_reaction_coordinate_projections(
        iguessmd_positions,
        TEST_iGUESSMD_LINEAR_PATH,
    )

    # Calculate reduced expected projections array (every nth point)
    reduced_expected_projections = expected_projections[:, ::every_nth]

    # Calculate reduced projections array (every nth point)
    reduced_calculated_projections = calculate_reaction_coordinate_projections(
        iguessmd_positions,
        TEST_iGUESSMD_LINEAR_PATH,
        every_nth_point=every_nth,
        include_end_point=False,
    )

    # Check results
    assert calculated_projections.shape == expected_projections.shape
    assert calculated_projections == pytest.approx(position_shifts[0])

    assert reduced_calculated_projections.shape == reduced_expected_projections.shape
    assert reduced_calculated_projections == pytest.approx(position_shifts[0])


@pytest.mark.parametrize(
    "reaction_path, every_nth",
    product(
        [TEST_iGUESSMD_LINEAR_PATH, TEST_iGUESSMD_CIRCULAR_PATH],
        [1, 2, 5, 23],
    )
)
def test_calculate_displacements_along_reaction_coordinate(reaction_path, every_nth):
    """
    Check that the displacements along the reaction path are
    correctly calculated for a given reaction path, including
    the cases where every_nth is greater than 1.
    """
    # Check results for full array
    expected_displacements = np.diff(reaction_path, axis=0)
    calculated_displacements = calculate_displacements_along_reaction_coordinate(reaction_path)
    assert (calculated_displacements == expected_displacements).all()

    # Check results for every_nth array
    reduced_expected_displacements = expected_displacements[::every_nth]
    reduced_calculated_displacements = calculate_displacements_along_reaction_coordinate(
        reaction_path,
        every_nth_point=every_nth,
        include_end_point=False
    )
    assert (reduced_calculated_displacements == reduced_expected_displacements).all()


@pytest.mark.parametrize(
    "reaction_path, every_nth",
    product(
        [TEST_iGUESSMD_LINEAR_PATH, TEST_iGUESSMD_CIRCULAR_PATH],
        [1, 2, 5, 23],
    )
)
def test_calculate_distance_along_reaction_coordinate(reaction_path, every_nth):
    """
    Check that the distance along the reaction coordinate is
    correctly calculated for a given reaction path, including
    the cases where every_nth is greater than 1.
    """
    # Check results for full array
    expected_distances = np.zeros(reaction_path.shape[0])
    expected_distances[1:] = np.cumsum(np.linalg.norm(np.diff(reaction_path, axis=0), axis=1), axis=0)
    calculated_distances = calculate_distance_along_reaction_coordinate(reaction_path)
    assert (calculated_distances == expected_distances).all()

    # Check results for every_nth array
    reduced_expected_distances = expected_distances[::every_nth]
    reduced_calculated_distances = calculate_distance_along_reaction_coordinate(
        reaction_path,
        every_nth_point=every_nth,
        include_end_point=False
    )
    assert (reduced_calculated_distances == reduced_expected_distances).all()


@pytest.mark.parametrize("n_steps", [100, 1000, 12345, 96235])
def test_calculate_pmf_exponential_average(n_steps):
    """
    Check that the PMF calculated for a set of example work arrays
    using the function calculate_pmf_exponential_average_kJ_mol is
    approximately equal to the expected PMF calculated explicitly.
    They may differ slightly numerically, as the
    """

    # Define temperature and calculate value of beta
    temp_K = 273.15
    beta_mol_kJ = calculate_beta_mol_kJ(temp_K)

    # Define example work arrays (array of zeros and monotonically
    # increasing arrays)
    work_array_1 = np.zeros(n_steps)
    work_array_2 = np.linspace(0,n_steps-1,n_steps)
    work_array_3 = np.linspace(0, 2*(n_steps - 1), n_steps)
    work_array_4 = np.linspace(0, 3*(n_steps - 1), n_steps)
    work_arrays = np.array([work_array_1, work_array_2, work_array_3, work_array_4])

    # Calculate expected exponential average PMF manually
    expected_pmf = - (1./beta_mol_kJ) * (np.log(np.sum(np.exp(-beta_mol_kJ * work_arrays), 0)) - np.log(work_arrays.shape[0]))

    # Calculate PMF via function
    calculated_pmf = calculate_pmf_exponential_average_kJ_mol(work_arrays, temp_K)

    # Check that expected and calculated values are close
    assert np.allclose(calculated_pmf, expected_pmf, rtol=1E-8)

    # Check that final value of calculated PMF is dominated by lowest work value (0),
    # which should yield approximately -(1 / beta) * log(1/N) = (1 / beta) * log(N)
    # [ONLY VALID FOR LOWER TEMPERATURES]
    assert calculated_pmf[-1] == (1./beta_mol_kJ) * np.log(work_arrays.shape[0])


@pytest.mark.parametrize(
    "mu, sigma",
    [
        (0.0, 1.0),
        (2.5, 0.3),
        (-0.5, 3.0),
    ],
)
def test_calculate_pmf_second_cumulant(mu, sigma):
    """
    Check that the PMF calculated for a set of Gaussian distributed
    work values using calculate_pmf_second_cumulant_kJ_mol yields
    the PMF expected for a Gaussian distribution.
    """
    # Define initial temperature and calculate beta
    temp_K = 273.15
    beta_mol_kJ = calculate_beta_mol_kJ(temp_K)

    # Define Gaussian distributed work arrays given mean mu
    # and standard deviation sigma
    n_samples = 10000000
    work_arrays = np.array([np.random.normal(mu, sigma, n_samples) for _ in range(10)]).transpose()

    # Vary means of distributions
    work_arrays += np.linspace(0, 9, 10)

    # Define expected PMF from second cumulant formula (Bessel-corrected variance)
    bessel_var = (sigma**2) * (n_samples)/(n_samples-1)
    expected_pmf = np.array([i + mu - (beta_mol_kJ/2.) * bessel_var for i in range(10)])

    # Calculate PMF via function
    calculated_pmf = calculate_pmf_second_cumulant_kJ_mol(work_arrays, temp_K)

    # Check that values are approximately equal (tolerance set
    # relatively high to account for sampling error)
    assert np.allclose(calculated_pmf, expected_pmf, rtol=5E-3)


def test_load_general_iguessmd_data():
    # TODO: Add test!
    pass


def test_load_iguessmd_simulation_data():
    # TODO: Add test!
    pass
