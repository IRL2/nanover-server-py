"""
Tests for :mod:`nanover.iguessmd.analysis`.

Things to test for analysis functions:
- General iGUESSMD data loads correctly [ ]
- Boltzmann constant is calculated correctly in units of kJ mol-1 K-1 [√]
- Beta (1 / (kB * T)) is correctly calculated in units of mol kJ-1 [√]
- Calculation of PMF via exponential average returns expected result [ ]
- Calculation of PMF via second cumulant approximation returns expected result [ ]
- Calculation of reaction coordinate projections works as expected [ ]
- Calculation of displacements along reaction path works as expected [ ]
- Calculation of distance along reaction coordinate works as expected [ ]
"""
import pytest
import numpy as np

from nanover.iguessmd.analysis import *

KB_KJ_MOL_K_VALUE = 0.008314462618

def test_calculate_boltzmann_constant_in_kJ_mol_K():
    """
    Check that the Boltzmann constant is correctly calculated in units of kJ mol-1 K-1
    """
    kB_calculated = boltzmann_constant_in_kJ_mol_K()
    assert kB_calculated._value == pytest.approx(KB_KJ_MOL_K_VALUE, rel=1E-10)


@pytest.mark.parametrize("temperature_K", [-1.0, 0.0, 1, np.pi, 273.15, 300.0, 12345.6789])
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
    beta_expected = 1. / (KB_KJ_MOL_K_VALUE * temperature_K)
    assert beta_calculated == pytest.approx(beta_expected, rel=1E-10)
