import itertools
from narupa.app import NarupaClient


def grouper(iterable, n, fillvalue=None):
    args = [iter(iterable)] * n
    return itertools.zip_longest(*args, fillvalue=fillvalue)


def atom_listing_from_frame(frame):
    neighbour_map = neighbour_map_from_frame(frame)

    atoms = []

    for index, element, id, residue in zip(itertools.count(0),
                                           frame.arrays['particle.element'],
                                           frame.arrays['atom.id'],
                                           frame.arrays['atom.residue']):
        atom = dict(index=index,
                    element=element,
                    id=id,
                    residue=residue,
                    neighbours=neighbour_map.get(index))
        atoms.append(atom)

    for atom in atoms:
        atom['neighbours'] = [atoms[index] for index in atom['neighbours']]

    return atoms


def neighbour_map_from_frame(frame):
    neighbour_map = {}

    for a, b in grouper(frame.arrays['bond'], 2):
        neighbours = neighbour_map.get(a, set())
        neighbours.add(b)
        neighbour_map[a] = neighbours

        neighbours = neighbour_map.get(b, set())
        neighbours.add(a)
        neighbour_map[b] = neighbours

    return neighbour_map


def main():
    with NarupaClient() as client:
        frame = client.wait_until_first_frame()
        atoms = atom_listing_from_frame(frame)

        selection = [atom['index'] for atom in atoms if atom['id'] == 'CG']
        selection = [atom for atom in atoms if len(atom['neighbours']) >= 3]
        print(selection)


if __name__ == '__main__':
    main()
