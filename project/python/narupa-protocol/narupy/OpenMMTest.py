import urllib.request
from time import sleep

from simtk.openmm.app import *

from narupa.protocol.topology.topology_pb2 import *
from narupy.StreamTopology import *


def download_pdb(pdb_id):
    url = 'http://www.rcsb.org/pdb/files/%s.pdb' % pdb_id
    pdb_text = urllib.request.urlretrieve(url, "%s.pdb" % pdb_id)
    return PDBFile("%s.pdb" % pdb_id)


# Prepare topology for sending to client
def convert_topology(topology):
    data = TopologyData()

    for residue in topology.residues():
        data.arrays['residue.id'].string_values.values.append(residue.name)
        data.arrays['residue.chain'].index_values.values.append(residue.chain.index)

    for atom in topology.atoms():
        data.arrays['atom.id'].string_values.values.append(atom.name)
        data.arrays['atom.element'].index_values.values.append(atom.element.atomic_number)
        data.arrays['atom.residue'].index_values.values.append(atom.residue.index)

    return data


from pdbfixer import PDBFixer

def load_and_fix_pdb(pdb_id, forcefield):
    print("Loading PDB")
    fixer = PDBFixer(filename="%s.pdb" % pdb_id)
    print("Find missing residues")
    fixer.findMissingResidues()
    print("Find non-standard residues")
    fixer.findNonstandardResidues()
    print("Replace non-standard residues")
    fixer.replaceNonstandardResidues()
    print("Remove heterogens")
    fixer.removeHeterogens(False)
    print("Find missing atoms")
    fixer.findMissingAtoms()
    print("Add missing atoms")
    fixer.addMissingAtoms()
    print("Add missing hydrogens")
    topology, positions = fixer.topology, fixer.positions
    modeller = Modeller(topology, positions)
    modeller.addHydrogens(forcefield, 7.0)
    return modeller.getTopology(), modeller.getPositions()


from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *


pdb = PDBFile("input.pdb")

topology, positions = pdb.topology, pdb.positions

server = NarupaServer()

client = NarupaClient()



forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')

system = forcefield.createSystem(topology, nonbondedMethod=PME,
        nonbondedCutoff=1*nanometer, constraints=HBonds)
integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
simulation = Simulation(topology, system, integrator)
simulation.context.setPositions(positions)
simulation.minimizeEnergy()
print("Step")
simulation.step(1)
