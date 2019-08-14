# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.
from typing import Optional, Collection

import numpy as np
import time
# from ase import Atoms
# from ase.calculators.calculator import Calculator, all_changes
from ctypes import CDLL, c_int, c_char_p, c_double, pointer
import ctypes
from numpy.ctypeslib import ndpointer

EV_PER_HARTREE = 27.2114
ANG_PER_BOHR = 0.529177
dynamics = 100

DFTBpluslib = CDLL('./libdftb+.so')

output_location = "tempfile.out"  # Output
output_mutable_string = ctypes.create_string_buffer(str.encode(output_location))

input_location = "dftb_in.hsd"  # Input
input_mutable_string = ctypes.create_string_buffer(str.encode(input_location))

func = DFTBpluslib.__getattr__("dftbp_init")
print(func)


def wrap_function(lib, funcname, restype, argtypes):
    """Simplify wrapping ctypes functions"""
    func = lib.__getattr__(funcname)
    func.restype = restype
    func.argtypes = argtypes
    return func


DFTB_state = ctypes.c_double(0.0)
DFTB_handler = ctypes.c_double(0.0)

wrap_function(DFTBpluslib, "dftbp_init", None, [ctypes.POINTER(ctypes.c_double), ctypes.c_char_p])
DFTBpluslib.dftbp_init(ctypes.pointer(DFTB_state), output_mutable_string)  # Initialize

wrap_function(DFTBpluslib, "dftbp_get_input_from_file", None, [ctypes.POINTER(ctypes.c_double), ctypes.c_char_p, ctypes.POINTER(ctypes.c_double)])
DFTBpluslib.dftbp_get_input_from_file(ctypes.pointer(DFTB_state), input_mutable_string, ctypes.pointer(DFTB_handler))  # Read input from file

wrap_function(DFTBpluslib, "dftbp_process_input", None, [ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double)])
DFTBpluslib.dftbp_process_input(ctypes.pointer(DFTB_state), ctypes.pointer(DFTB_handler))  # Process input

wrap_function(DFTBpluslib, "dftbp_get_nr_atoms", ctypes.c_int, [ctypes.POINTER(ctypes.c_double)])
natom = DFTBpluslib.dftbp_get_nr_atoms(ctypes.pointer(DFTB_state))  # Get number of atoms
print(natom)
print('\n')

# Generate a general numpy matrix that is c compatible
_doublepp = ndpointer(dtype=np.float64, ndim=1, flags='C')
x = np.arange(natom * 3)
xyz = (x.__array_interface__['data'][0] + np.arange(x.shape[0]) * x.strides[0]).astype(np.float64)
forces = (x.__array_interface__['data'][0] + np.arange(x.shape[0]) * x.strides[0]).astype(np.float64)

start_time = time.time()
for i in range(0, dynamics):

    #print(i)
    #print('\n')

    wrap_function(DFTBpluslib, "dftbp_get_coords", None, [ctypes.POINTER(ctypes.c_double), _doublepp])
    DFTBpluslib.dftbp_get_coords(ctypes.pointer(DFTB_state), xyz)  # Get coordinates
    for i in range(0, len(xyz)):
        xyz[i] = xyz[i]/1.889725989
    #print(xyz)
    #print('\n')

    wrap_function(DFTBpluslib, "dftbp_get_gradients", None, [ctypes.POINTER(ctypes.c_double), _doublepp])
    DFTBpluslib.dftbp_get_gradients(ctypes.pointer(DFTB_state), forces)  # Get forces
    #print(forces)
    #print('\n')

    for i in range(0, len(xyz)):
        xyz[i] = xyz[i]*1.889725989 + 0.01  # Modify the coordinates

    wrap_function(DFTBpluslib, "dftbp_set_coords", None, [ctypes.POINTER(ctypes.c_double), _doublepp])
    DFTBpluslib.dftbp_set_coords(ctypes.pointer(DFTB_state), xyz)  # Set coordinates
    #print("managed to load without crashing")

print("Time = --- %s seconds ---" % (time.time() - start_time))

# class DFTBCalculator(Calculator):
#     """
#     Simple implementation of an ASE calculator for Sparrow.
#
#     Parameters:
#         method :  The electronic structure method to use in calculations.
#     """
#     implemented_properties = ['energy', 'forces']
#
#     def __init__(self, atoms: Optional[Atoms] = None, method='DFTB0', **kwargs):
#         super().__init__(**kwargs)
#         self.atoms = atoms
#         self.method = method
#         self.DFTBpluslib = CDLL('./libdftb+.so')
#
#
#         output_location = c_char_p("tempfile.out")
#         DFTB_pointer = pointer(c_double)
#         self.DFTBpluslib.dftbp_init(DFTB_pointer, output_location)
#
#
#     def calculate(self, atoms: Optional[Atoms] = None,
#                   properties=('energy', 'forces'),
#                   system_changes=all_changes):
#
#         if atoms is None:
#             atoms = self.atoms
#             raise ValueError('No ASE atoms supplied to calculator, and no ASE atoms supplied with initialisation.')
#         self._calculate_sparrow(atoms, properties)
#
#     def _calculate_sparrow(self, atoms: Atoms, properties: Collection[str]):
#         positions = atoms.positions
#         elements = atoms.get_chemical_symbols()
#
#         kwargs = {property_name: True for property_name in properties}
#         # TODO pass these to calculate in wrapper.
#         if 'energy' in properties:
#             energy_hartree = calculate_energy(elements, positions, self.method)
#             self.results['energy'] = energy_hartree * EV_PER_HARTREE
#         if 'forces' in properties:
#             #TODO make np array come out of wwrapper.
#             gradients_hartree_bohr = np.array(calculate_gradients(elements, positions, self.method))
#             self.results['forces'] = - gradients_hartree_bohr * EV_PER_HARTREE / ANG_PER_BOHR
#         return


