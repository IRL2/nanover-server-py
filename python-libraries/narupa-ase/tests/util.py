from ase import Atoms


def co_atoms():
    d = 1.1
    co = Atoms('CO', positions=[(0, 0, 0), (0, 0, d)],
               cell=[20, 20, 20],
               pbc=[1, 1, 1])
    return co
