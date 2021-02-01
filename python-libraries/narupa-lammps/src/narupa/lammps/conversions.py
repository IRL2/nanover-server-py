# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.
"""
Module providing conversion utilities between Narupa and LAMMPS.
"""

from typing import NamedTuple

# LAMMPS works with arbitrary masses, so we need to convert it to a nuclear number
# This list is a best guess for atom types, but won't work for isotopes for now.
# see issue 82 https://gitlab.com/intangiblerealities/narupa-protocol/issues/82
ELEMENT_INDEX_MASS = {
    1: 1,
    3: 1,
    4: 2,
    7: 3,
    9: 4,
    11: 5,
    12: 6,
    14: 7,
    16: 8,
    19: 9,
    20: 10,
    23: 11,
    24: 12,
    27: 13,
    28: 14,
    31: 15,
    32: 16,
    35: 17,
    39: 19,
    40: 18,
    45: 21,
    48: 22,
    51: 23,
    52: 24,
    55: 25,
    56: 26,
    59: 27,
    64: 29,
    65: 30,
    70: 31,
    73: 32,
    75: 33,
    79: 34,
    80: 35,
    84: 36,
    85: 37,
    88: 38,
    89: 39,
}

# Check what units are being used in LAMMPS using this dict
# For now support converting lj and real, for full unit list
# see https://lammps.sandia.gov/doc/units.html
# List consists of the unit type, the conversion required to convert positions to nm,
# and the conversion required to convert forces to kJ/mol/nm
class LammpsUnitConverter(NamedTuple):
    type: str
    positions: float
    forces: float

# store plank values as a list so that we don't do float lookups in a dict.
LAMMPS_UNITS_CHECK = {
    # Lennard jones: Is unitless, everything is set to 1
    0: LammpsUnitConverter(type="lj", positions=1, forces=1),
    # Real:
    # Distance: 1 angstrom- > nm (10)
    # Force:    kj/mol/angstrom -> kcal/mol/nm (4.1840 *10) (Confirmed in MDanaysis)
    1: LammpsUnitConverter(type="real", positions=10, forces=41.840),
    # Metal:
    # Distance: angstrom -> nm, (10)
    # Force: eV/angstrom -> kcal/mol/nm (96.485*10) (Confirmed in MDanalysis)
    2: LammpsUnitConverter(type="metal", positions=10, forces=964.85),
    # SI:
    # Distance: meters ->nm (10^-9)
    # Force: Newtons q-> kcal/mol/nm (602214128999.9999)
    3: LammpsUnitConverter(type="si", positions=10 ** -9, forces=602214128999.9999),
    # cgs:
    # Distance: centemters -> nm
    # Froce: dyne (1/100000 Newtons) -> kj/mol*nm
    4: LammpsUnitConverter(type="cgs", positions=10 ** -7, forces=6022141.289999999),
    # Electron:
    # Distance: Bohr -> nm
    # Force: Hartree/Bohr (2625.50 / 0.0529117) ->  kj/mol*nm
    5: LammpsUnitConverter(type="electron", positions=0.05529177, forces=49620.4053),
    # Mirco:
    # Distance: mircometers -> nm,
    # Force: pircogram-micrometer/microsecond^2 -> Newtons
    # (1/1000000000000 *((1/1000000)/(1/1000000)^2) ->  kj/mol*nm
    6: LammpsUnitConverter(type="micro", positions=1000, forces=60221.41289999999),
    # Nano:
    # Distance: nanometers,
    # Force: atoogram-nanometer/nanosecond^2  -> Newtons
    # (1/1e-12 *((1/1e-9)/(1/1e-9)^2) ->  kj/mol*nm
    7: LammpsUnitConverter(type="nano", positions=1.0, forces=602214128.9999999)
}
PLANK_VALUES = (
    1.0,
    95.306976368,
    4.135667403e-3,
    6.62606896e-34,
    6.62606896e-27,
    0.1519829846,
    6.62606896e-13,
    6.62606896e-4
)