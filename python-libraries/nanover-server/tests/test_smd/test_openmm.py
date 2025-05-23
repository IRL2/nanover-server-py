"""
Tests for :mod:`nanover.smd.openmm`.

Things to test:
- An SMD simulation can be created either from an existing OpenMM simulation or
  a NanoVer OpenMM XML file
- OpenMMSMDSimulation returns OpenMMSMDSimulationAtom or OpenMMSMDSimulationCOM
  as appropriate
- The PBCs of the loaded simulation are respected by the SMD force, and the SMD
  force shares this periodicity
- The SMD force attaches to the correct atom (dictated by the index/indices passed
  to the class upon creation)
- The simulation can be reset correctly, with all attributes returning to the same
  state as immediately after the creation of the class itself
- SMD force is correctly added to the system
- SMD force can be correctly removed from the system
- SMD force position is correctly updated
- Running the equilibration with the initial restraint throws an error correctly
  when run with the SMD force not located at the initial position
- The SMD simulation can be saved correctly, with or without the SMD force
- If the SMD simulation is being loaded from a NanoVer OpenMM XML file created
  via one of the SMDSimulation classes and contains an SMD force already, this
  SMD force is correctly loaded and matches the expected force constant specified
  when creating the class
- The class can generate the correct number of starting structures in the specified
  time interval, and that these are saved to the correct location
- Running an SMD simulation produces reasonable results for the cumulative work done
  (may need to think about a specific test case for this...)
- _calculate_forces works as expected
- _calculate_work_done works as expected
- Simulation data is saved in the correct format to the correct location, and can be
  subsequently loaded back into python correctly
- General SMD data is saved in the correct format to the correct location, and can be
  subsequently loaded back into python correctly
- For OpenMMSMDSimulationCOM, the COM of the specified atoms is correctly calculated
- For OpenMMSMDSimulationCOM, the trajectory of the COM of the specified atoms is
  correctly calculated
- smd_com_force works as expected
- smd_single_atom_force works as expected
"""

import openmm as mm

from nanover.smd.openmm import *

from ..test_openmm.simulation_utils import (
    basic_system,
    basic_simulation,
    basic_simulation_with_imd_force,
    BASIC_SIMULATION_POSITIONS,
    empty_imd_force,
    assert_basic_simulation_topology,
    single_atom_system,
    single_atom_simulation,
    single_atom_simulation_with_imd_force,
    ARGON_SIMULATION_POSITION,
    assert_single_atom_simulation_topology,
)

# Very basic thing to test entire class as it would be used: tutorial notebook that can be tested