from ase import Atoms
import itertools
import numpy as np

from narupa.protocol.trajectory import FrameData


AngToNm = 0.1
NmToAng = 1.0 / AngToNm

atom_radiuses_ang = {"H": 1.2, "Cu": 2.38, "Pt": 2.32}

def ase_to_framedata(ase_atoms: Atoms, positions=True, topology=True, state=True) -> FrameData:
    data = FrameData()

    if topology:
        data.arrays['residue.id'].string_values.values.extend(["ASE"])
        data.arrays['residue.chain'].index_values.values.extend([0])

        for index, atom in enumerate(ase_atoms):
            data.arrays['atom.id'].string_values.values.append(str(atom.index))
            data.arrays['atom.element'].index_values.values.append(atom.number)
            data.arrays['atom.residue'].index_values.values.append(0)

        bonds = GenerateBonds(ase_atoms)
        for bond in bonds:
            data.arrays['bond'].index_values.values.append(bond[0])
            data.arrays['bond'].index_values.values.append(bond[1])

    if positions:
        array = data.arrays['atom.position'].float_values.values
        floats = ase_atoms.get_positions().flatten() * AngToNm
        array[:] = floats

    if state:
        data.values['energy.potential'].number_value = ase_atoms.get_potential_energy()
        data.values['energy.kinetic'].number_value = ase_atoms.get_kinetic_energy()

    return data


def GetRadius(symbol):
    """
    Gets the radius of an atom (Angstrom)
    :param symbol:
    :return:
    """
    return atom_radiuses_ang[symbol]


def Threshold(radii):
    return 0.6 * sum(radii)


def GenerateBonds(atoms: Atoms):
    bonds = []
    for pair in itertools.combinations(atoms, 2):
        distance = np.linalg.norm(pair[0].position - pair[1].position)
        symbol = pair[0].symbol
        radii = [GetRadius(atom.symbol) for atom in pair]
        if distance < Threshold(radii):
            bonds.append([atom.index for atom in pair])
    return bonds
