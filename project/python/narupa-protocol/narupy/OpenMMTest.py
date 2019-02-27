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

    for residue in topology.residues():
        data.arrays['residue.id'].string_values.values.append(residue.name)
        data.arrays['residue.chain'].index_values.values.append(residue.chain.index)

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

pdbs = ['1A9X','1AON','1AUY','1B35','1BEV','1BR2','1BXR','1C2W','1C30','1C3O','1CE8','1CS0','1CWP','1CX8','1D4M','1DEQ','1DWN','1DZL','1E7P','1E9S','1EI7','1EJ6','1F1H','1F2N','1F52','1F8V','1FFK','1FJG','1FKA','1FNT','1FPN','1FPY','1G0U','1G3I','1G65','1GAV','1GQ2','1GR5','1GRU','1GT7','1GW7','1GW8','1GYT','1H1K','1H2I','1H6D','1H8T','1HB7','1HB9','1HNW','1HNX','1HNZ','1HR0','1HTO','1HTQ','1I41','1I43','1I48','1I94','1I95','1I96','1I97','1IBK','1IBL','1IBM','1IR2','1IRU','1IW7','1IWA','1IZL','1J0B','1J5A','1J5E','1JD2','1JDB','1JGO','1JGP','1JGQ','1JJ2','1JRO','1JRP','1JS9','1JZX','1JZY','1JZZ','1K01','1K32','1K3V','1K5M','1K73','1K8A','1K9M','1KC8','1KD1','1KEE','1KP8','1KQS','1KYI','1KYO','1L9U','1LAJ','1LGR','1LP3','1M1C','1M1K','1M1Y','1M34','1M6V','1M8Q','1M90','1MCZ','1MJG','1ML5','1MNF','1MQT','1MT5','1MVW','1MX9','1N32','1N33','1N34','1N36','1N4P','1N4Q','1N4R','1N4S','1N6D','1N6E','1N6F','1N8R','1NG0','1NJI','1NJM','1NJN','1NJO','1NJP','1NKW','1NQT','1NR7','1NT9','1NWX','1NWY','1O18','1O19','1O1A','1O1B','1O1C','1O1D','1O1E','1O1F','1O1G','1OGY','1OHF','1OHG','1OJD','1OND','1OPO','1P5W','1P5Y','1P7G','1P9X','1PCQ','1PDI','1PDP','1PF9','1PGL','1PGW','1PIV','1PJL','1PMA','1PQV','1PVC','1Q3G','1Q7Y','1Q81','1Q82','1Q86','1QGT','1QJU','1QJX','1QJY','1QO5','1QOX','1QVF','1QVG','1QZV','1RBL','1RCO','1RCX','1RSC','1RVV','1RXS','1RYP','1RYY','1S5L','1S64','1S72','1SID','1SIE','1SM1','1SMY','1SOJ','1SVA','1SVT','1SX3','1SX4','1T36','1T5E','1T60','1T6P','1TNB','1TNO','1TNU','1TNY','1TNZ','1TZL','1TZN','1TZO','1U1Y','1UF2','1UPM','1UW9','1UWA','1UZD','1UZH','1V9U','1VAK','1VB2','1VB4','1VBA','1VBB','1VBC','1VBE','1VF7','1VPX','1VQ4','1VQ5','1VQ6','1VQ7','1VQ8','1VQ9','1VQK','1VQL','1VQM','1VQN','1VQO','1VQP','1VVJ','1VY4','1VY5','1VY6','1VY7','1VYH','1W1I','1W2B','1W36','1W4C','1W5C','1W63','1W8X','1WCD','1WCE','1WCM','1WE5','1X33','1X35','1X36','1X9P','1X9T','1XBP','1XCK','1XFU','1XFV','1XFW','1XFX','1XFY','1XFZ','1XI4','1XI5','1XMO','1XMQ','1XNQ','1XNR','1XSI','1XSJ','1XSK','1Y0V','1Y1V','1Y1W','1Y1Y','1Y69','1Y77','1YA7','1YAJ','1YAR','1YAU','1YG8','1YHQ','1YI2','1YIJ','1YIT','1YJ9','1YJN','1YJW','1YNJ','1YNN','1YQ2','1YRQ','1Z14','1Z58','1Z6O','1Z7Q','1ZA7','1ZBA','1ZBE','1ZKU','1ZUM','1ZY8','1ZYR','2A68','2A69','2A6E','2A6H','2AAR','2AAZ','2AFI','2AXT','2B2D','2B2E','2B2G','2B63','2B8K','2B9V','2BBV','2BE5','2BGN','2BNY','2BQ5','2BR2','2BS1','2BT4','2BTV','2BU1','2BUK','2C37','2C38','2C39','2C4Q','2C4Y','2C4Z','2C50','2C51','2C6S','2C7C','2C7D','2C7E','2CB6','2CDH','2CGT','2CKJ','2CSE','2CW0','2D3O','2DF7','2DW7','2E0Z','2E1Q','2E5L','2EU1','2EWO','2EX3','2F16','2F2H','2F4V','2F5Z','2FAK','2FHG','2FHH','2FJF','2FK0','2FL8','2FL9','2FRV','2FTC','2FUG','2FYN','2FZ1','2FZ2','2G33','2G34','2G8G','2GH8','2GLJ','2GLS','2GPL','2GRE','2GSY','2GTL','2GTT','2H1L','2HG4','2HHH','2HLD','2HRT','2I2X','2IGK','2IGM','2IGN','2IGO','2IJZ','2IZ8','2IZ9','2IZN','2IZW','2J28','2J7A','2JA5','2JA6','2JA7','2JA8','2JD6','2JD7','2JD8','2JES','2JIZ','2JJ1','2JJ2','2JJM','2LGS','2NOX','2NV2','2NWC','2O01','2O1T','2O43','2O44','2O45','2O5I','2O5J','2OGM','2OGN','2OGO','2ONM','2OTJ','2OTL','2PFF','2PMZ','2PNK','2PPB','2Q08','2Q3E','2QA0']
#download_and_fix_pdb(pdb_id, forcefield)

for pdb_id in pdbs[100:]:

    try:
        pdb = download_pdb(pdb_id)

        topology, positions = pdb.topology, pdb.positions

        if(len(positions) < 30000):
            continue

        for i in range(1, 100):
            convert_frame(positions)

        timings = []

        for i in range(1, 1000):
            start_time = time.time()
            convert_frame(positions)
            end_time = time.time()
            timings.append(end_time - start_time)

        import statistics

        print("%s, %i, %f, %f, %f" % (pdb_id, len(positions), statistics.median(timings), statistics.mean(timings), statistics.stdev(timings)))
    except:
        pass
exit(0)

server = NarupaServer()

#client = NarupaClient()

def on_frame(list):
    print("FRAME")

def on_topology(list):
    print("TOPOLOGY")


system = forcefield.createSystem(topology, nonbondedMethod=PME,
        nonbondedCutoff=1*nanometer, constraints=HBonds)
integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
simulation = Simulation(topology, system, integrator)
simulation.context.setPositions(positions)
simulation.minimizeEnergy()


i = 0

server.send_topology(i, convert_topology(simulation.topology))

start_time = time.time()

while True:
    simulation.step(1)
    convert_frame(simulation.context.getState(getPositions=True).getPositions())
    #server.send_frame(i, convert_frame(simulation.context.getState(getPositions=True).getPositions()))
    i += 1
    if i % 60 == 0:
        print(i / (time.time() - start_time))
