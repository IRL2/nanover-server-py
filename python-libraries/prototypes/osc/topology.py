# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.
"""
Functions for extracting topology from narupa frames into per-atom dictionaries.
"""
from narupa.trajectory.frame_data import PARTICLE_ELEMENTS, PARTICLE_RESIDUES, PARTICLE_NAMES, BOND_PAIRS

ATOM_PROPERTIES = {
    'element': PARTICLE_ELEMENTS,
    'name': PARTICLE_NAMES,
    'residue': PARTICLE_RESIDUES,
}


def atom_listing_from_frame(frame):
    """
    Create a list with one dictionary per atom with entries, where present, for
    the id, element, residue, and neighbour information of each atom.
    """
    atom_arrays = {}

    for name, key in ATOM_PROPERTIES.items():
        array = frame.arrays.get(key)
        if array is not None:
            atom_arrays[name] = array

    if not atom_arrays:
        return []

    atom_count = min(len(array) for array in atom_arrays.values())
    atoms = []
    neighbour_map = neighbour_map_from_frame(frame)

    for index in range(0, atom_count):
        atom = {'index': index, 'neighbours': neighbour_map.get(index)}

        for key, array in atom_arrays.items():
            atom[key] = array[index]

        atoms.append(atom)

    for atom in atoms:
        if atom['neighbours'] is not None:
            atom['neighbours'] = [atoms[index] for index in atom['neighbours']]

    return atoms


def neighbour_map_from_frame(frame):
    """
    Map out neighbour relationships of atoms from the bond information present
    in a narupa frame.
    """
    neighbour_map = {}

    if BOND_PAIRS not in frame.arrays:
        return neighbour_map

    for a, b in frame.bond_pairs:
        neighbours = neighbour_map.get(a, set())
        neighbours.add(b)
        neighbour_map[a] = neighbours

        neighbours = neighbour_map.get(b, set())
        neighbours.add(a)
        neighbour_map[b] = neighbours

    return neighbour_map
