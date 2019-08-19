"""
Functions for extracting topology from narupa frames into per-atom dictionaries.
"""

ATOM_PROPERTIES = {
    'element': 'particle.element',
    'id': 'atom.id',
    'residue': 'atom.residue',
}


def grouper(iterable, n):
    args = [iter(iterable)] * n
    return zip(*args)


def atom_listing_from_frame(frame):
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
    neighbour_map = {}

    if 'bond' not in frame.arrays:
        return neighbour_map

    for a, b in grouper(frame.arrays['bond'], 2):
        neighbours = neighbour_map.get(a, set())
        neighbours.add(b)
        neighbour_map[a] = neighbours

        neighbours = neighbour_map.get(b, set())
        neighbours.add(a)
        neighbour_map[b] = neighbours

    return neighbour_map
