import time
from sys import stdout

from openmm.app import *
from openmm import *
from openmm.unit import *

from nanover.imd import ParticleInteraction
from nanover.omni import OmniRunner
from nanover.omni.openmm import OpenMMSimulation
from test_dcd import dcd_reporter

psf = CharmmPsfFile('smd-test/da.psf')
pdb = PDBFile('smd-test/smd_ini.pdb')
params = CharmmParameterSet('smd-test/top_all27_prot_lipid.rtf', 'smd-test/par_all27_prot_lipid.prm')

system = psf.createSystem(params, nonbondedMethod=NoCutoff, constraints=HBonds,  hydrogenMass=1.5*amu)

restraint = CustomExternalForce('k*periodicdistance(x, y, z, x0, y0, z0)^2')
system.addForce(restraint)
restraint.addGlobalParameter('k', 3012.48*kilojoules_per_mole/nanometer**2)
restraint.addPerParticleParameter('x0')
restraint.addPerParticleParameter('y0')
restraint.addPerParticleParameter('z0')
restraint.addParticle(1, pdb.positions[1]) # add constrain to N

integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.004*picoseconds)
simulation = Simulation(psf.topology, system, integrator)
simulation.context.setPositions(pdb.positions)
simulation.minimizeEnergy()

try:
    simulation.step(1000)

    dcd_reporter = DCDReporter('output1.dcd', 10)
    simulation.reporters.append(dcd_reporter)
    simulation.reporters.append(StateDataReporter(stdout, 200, step=True, potentialEnergy=True, temperature=True))

    omm = OpenMMSimulation.from_simulation(simulation)
    omm.frame_interval = 50

    with OmniRunner.with_basic_server(omm) as runner:
        runner.load(0)
        runner.app_server.imd.insert_interaction("interaction.test", ParticleInteraction(
            particles=[len(pdb.positions) - 4],
            position=[0, 0, 10],
            interaction_type="constant",
            scale=16,
        ))

        while True:
            time.sleep(.1)
finally:
    dcd_reporter.__del__()
