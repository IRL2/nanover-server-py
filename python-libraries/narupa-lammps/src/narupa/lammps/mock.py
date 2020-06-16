# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.
"""
Module providing a mock LAMMPS object so that the Narupa LAMMPS flow can be tested
without LAMMPS installed.
"""
import ctypes
from typing import List, Union


class MockLammps:
    """
    A fake LAMMPS object intended just for unit testing the LAMMPS code
    without having to have LAMMPS installed on a server
    """

    def __init__(self, n_atoms_in_dummy: int = None):
        # Set a default atom length for tests
        _DEFAULT_ATOMS = 3
        self.n_atoms = n_atoms_in_dummy if n_atoms_in_dummy is not None else _DEFAULT_ATOMS

    def gather_atoms(self, array_type: str, _dummy_variable, _array_shape):
        """
        Generates fake ctypes to mimic lammps internal pointers

        :param array_type: determines the type of data that should be replicated
        :param _dummy_variable: Unused here, only relevant to lammps
        :param _array_shape: Unused here, only relevant to lammps
        :return: matrix data_array that contains all the mock data
        """
        data_array: Union[ctypes.Array[ctypes.c_double], ctypes.Array[ctypes.c_int]]
        empty_list: List = []
        if array_type == "x":
            data_array = (ctypes.c_double * (3 * self.n_atoms))(*range(3 * self.n_atoms))
        elif array_type == "f":
            data_array = (ctypes.c_double * (3 * self.n_atoms))(*empty_list)
        elif array_type == "type":
            data_array = (ctypes.c_int * self.n_atoms)(*empty_list)
            # All atoms have the same type for DummyLammps testing
            for i in range(self.n_atoms):
                data_array[i] = 1
        else:
            raise Exception('Unknown array type asked for in dummyLammps.gather_atoms')

        return data_array

    def scatter_atoms(self, _array_type, _dummy_variable, _array_shape, __data_array):
        """
        Mimics lammp_class.scatter_atoms, in the mock case it does nothing
        Note: This can't be static as it is designed to replicate the internal LAMMPS class
        """
        return None

    def extract_global(self, types: str, _number_type):
        """
        Generate mock element list for testing

        replicates lammp_class.extract_global("ntypes", 0)
        :param types: LAMMPS global variable that is needed
        :param _number_type: unused
        :return: Element list.
        """
        dummy_element_list: Union[float, int]

        if types == "ntypes":
            dummy_element_list = 1
        elif types == "hplanck":
            dummy_element_list = 95.306976368
        else:
            raise Exception('Unknown array type asked for in dummyLammps.extract_global')

        return dummy_element_list

    def get_natoms(self):
        """
        Generate mock element list for testing
        """
        return self.n_atoms

    def extract_atom(self, types: str, _number_type):
        """
        Generate mock element list for testing
        :param types: string that indicates the info that should be passed
        :param _number_type: unused parameter to indicate integer float. etc
        :return: dummy_element_list
        """
        if types == "mass":
            # initialise mock type
            dummy_element_type = ctypes.c_double * 2
            # great array of mock type
            dummy_element_list = dummy_element_type()
            # For some reason masses have a blank [0] value in LAMMPS
            dummy_element_list[0] = 0
            dummy_element_list[1] = 1
        else:
            raise Exception('Unknown array type asked for in dummyLammps.extract_atom')

        return dummy_element_list
