import urllib.request
from time import sleep

from simtk.openmm.app import *

from narupa.protocol.topology.topology_pb2 import *
from narupy.StreamTopology import *


def download_pdb(pdb_id):
    url = 'http://www.rcsb.org/pdb/files/%s.pdb' % pdb_id
    urllib.request.urlretrieve(url, "%s.pdb" % pdb_id)
    return PDBFile("%s.pdb" % pdb_id)


# Prepare topology for sending to client
def convert_topology(topology):
    data = TopologyData()

    data.arrays['residue.id'].string_values.values.extend([residue.name for residue in topology.residues()])
    data.arrays['residue.chain'].index_values.values.extend([residue.chain.index for residue in topology.residues()])

    for atom in topology.atoms():
        data.arrays['atom.id'].string_values.values.append(atom.name)
        data.arrays['atom.element'].index_values.values.append(atom.element.atomic_number)
        data.arrays['atom.residue'].index_values.values.append(atom.residue.index)

    for bond in topology.bonds():
        data.arrays['bond'].index_values.values.append(bond[0].index)
        data.arrays['bond'].index_values.values.append(bond[1].index)

    return data

# Prepare topology for sending to client
def convert_frame(positions):
    data = FrameData()
    array = data.arrays['atom.position'].float_values.values
    floats = [value for position in positions for value in position._value]
    array.extend(floats)
    return data

from pdbfixer import PDBFixer

def download_and_fix_pdb(pdb_id, forcefield):
    url = 'http://www.rcsb.org/pdb/files/%s.pdb' % pdb_id
    urllib.request.urlretrieve(url, "%s.pdb" % pdb_id)
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
    PDBFile.writeFile(modeller.getTopology(), modeller.getPositions(), open("%s.pdb" % pdb_id, "w+"))


from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *

forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')

def test_timings():
    pdbs = ['5ZJZ','5ZK5','5ZK7','5ZK9','5ZML','5ZNL','5ZNS','5ZNT','5ZZE','6A10','6A50','6A53','6A54','6A55','6A6M','6A8G','6A8N','6AKF','6AKG','6BQD','6CCI','6CDT','6CPT','6CPU','6CT8','6CU2','6CXB','6CXG']

    for pdb_id in pdbs:

        #try:
            pdb = download_pdb(pdb_id)

            topology, positions = pdb.topology, pdb.positions

            for i in range(1, 100):
                convert_frame(positions)

            timings = []

            for i in range(1, 100):
                start_time = time.time()
                for i in range(1, 100):
                    convert_frame(positions)
                end_time = time.time()
                timings.append((end_time - start_time)/100)

            data = convert_frame(positions)

            import statistics

            print("%s, %i, %i, %f, %f, %f" % (pdb_id, len(positions), len(data.SerializeToString()), statistics.median(timings), statistics.mean(timings), statistics.stdev(timings)))
        #except:
         #   pass
    exit(0)

download_and_fix_pdb("5ZJZ", forcefield)
pdb = PDBFile("5ZJZ.pdb")
topology, positions = pdb.topology, pdb.positions

server = NarupaServer()

#client = NarupaClient()

def on_frame(list):
    print("FRAME")

def on_topology(list):
    print("TOPOLOGY")

system = forcefield.createSystem(topology, nonbondedMethod=PME,
        nonbondedCutoff=1*nanometer, constraints=HBonds)
system.AddInteractiveForce()
integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
simulation = Simulation(topology, system, integrator)
simulation.context.setPositions(positions)
simulation.minimizeEnergy()

i = 0

server.send_topology(i, convert_topology(simulation.topology))

start_time = time.time()

print("Starting simulation")
while True:
    simulation.step(1)
    convert_frame(simulation.context.getState(getPositions=True).getPositions())
    server.send_frame(i, convert_frame(simulation.context.getState(getPositions=True).getPositions()))
    i += 1
