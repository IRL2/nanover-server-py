from simtk.openmm.app import *

forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')

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