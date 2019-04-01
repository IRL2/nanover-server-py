from narupa.protocol.topology.topology_pb2 import TopologyData
from ase import Atoms
import itertools
import numpy as np

atom_radiuses_ang = {"H": 1.2, "Cu" : 2.38, "Pt" :  2.32 }
AngToNm = 0.1

def ase_atoms_to_topology_data(ase_atoms: Atoms) -> TopologyData:
    data = TopologyData()

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


