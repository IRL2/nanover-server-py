from ase import Atoms
import itertools
import numpy as np

from narupa.trajectory import FrameData

AngToNm = 0.1


def ase_atoms_to_frame_data(ase_atoms: Atoms) -> FrameData:
    data = FrameData()

    return data


atom_radiuses_ang = {"H": 1.2, "Cu": 2.38, "Pt": 2.32}
AngToNm = 0.1


def ase_to_framedata(ase_atoms: Atoms, positions=True, topology=True, state=True) -> FrameData:
    data = FrameData()

    if topology:
        data.arrays['residue.id'] = ["ASE"]
        data.arrays['residue.chain'] = [0]

        atom_names = []
        elements = []
        residue_ids = []
        for index, atom in enumerate(ase_atoms):
            atom_names.append(str(atom.index))
            elements.append(atom.number)
            residue_ids.append(0)

        data.arrays['atom.id'] = atom_names
        data.particle_elements = elements
        data.arrays['atom.residue'] = residue_ids

        bonds = GenerateBonds(ase_atoms)
        data.bonds = bonds

    if positions:
        data.particle_positions = ase_atoms.get_positions() * AngToNm

    if state:
        data.values['energy.potential'] = ase_atoms.get_potential_energy()
        data.values['energy.kinetic'] = ase_atoms.get_kinetic_energy()

    return data


def GetRadius(symbol):
    """
    Gets the radius of an atom (Angstrom)
    :param symbol:
    :return:
    """
    return atom_radiuses_ang[symbol]


def Threshold(radiuses):
    return 0.6 * sum(radiuses)


def GenerateBonds(atoms: Atoms):
    bonds = []
    for pair in itertools.combinations(atoms, 2):
        distance = np.linalg.norm(pair[0].position - pair[1].position)
        symbol = pair[0].symbol
        radiuses = [GetRadius(atom.symbol) for atom in pair]
        if distance < Threshold(radiuses):
            bonds.append([atom.index for atom in pair])
    return bonds
