import pytest

import numpy as np

import simtk.openmm as mm
from simtk.openmm import app
from simtk.unit import kelvin, picosecond, femtosecond, nanometer, amu


class DoNothingReporter:
    """
    OpenMM reporter that does nothing.

    The reporter does nothing but is valid. It is meant to populate the list of
    reporters of an OpenMM simulation.
    """
    # The name of the method is part of the OpenMM API. It cannot be made to
    # conform PEP8.
    def describeNextReport(self, simulation):  # pylint: disable=invalid-name
        """
        Activate the reporting every step, but collect no data.
        """
        return 1, False, False, False, False

    def report(self, simulation, state):
        """
        Do not report anything.
        """
        pass


@pytest.fixture
def basic_simulation():
    """
    Setup a minimal OpenMM simulation with two methane molecules.
    """
    periodic_box_vector = [
        [50,  0,  0],
        [ 0, 50,  0],
        [ 0,  0, 50]
    ]
    positions = np.array([
        # First residue
        [ 0,       0,      0],  # C
        [ 5.288,   1.610,  9.359],  # H
        [ 2.051,   8.240, -6.786],  # H
        [-10.685, -0.537,  1.921],  # H
        # Second residue, copied from the first but shifted
        # by 5 nm along the Z axis
        [  0,      0,      5],  # C
        [  5.288,  1.610, 14.359],  # H
        [  2.051,  8.240, -1.786],  # H
        [-10.685, -0.537,  6.921],  # H
    ], dtype=np.float32)

    topology = app.Topology()
    carbon = app.Element.getBySymbol('C')
    hydrogen = app.Element.getBySymbol('H')
    chain = topology.addChain()
    residue = topology.addResidue(name='METH1', chain=chain)
    c1 = topology.addAtom(element=carbon, name='C1', residue=residue)
    h2 = topology.addAtom(element=hydrogen, name='H2', residue=residue)
    h3 = topology.addAtom(element=hydrogen, name='H3', residue=residue)
    h4 = topology.addAtom(element=hydrogen, name='H4', residue=residue)
    topology.addBond(c1, h2)
    topology.addBond(c1, h3)
    topology.addBond(c1, h4)
    chain = topology.addChain()
    residue = topology.addResidue(name='METH2', chain=chain)
    c1 = topology.addAtom(element=carbon, name='C1', residue=residue)
    h2 = topology.addAtom(element=hydrogen, name='H2', residue=residue)
    h3 = topology.addAtom(element=hydrogen, name='H3', residue=residue)
    h4 = topology.addAtom(element=hydrogen, name='H4', residue=residue)
    topology.addBond(c1, h2)
    topology.addBond(c1, h3)
    topology.addBond(c1, h4)

    system = mm.System()
    system.setDefaultPeriodicBoxVectors(*periodic_box_vector)
    system.addParticle(mass=12)
    system.addParticle(mass=1)
    system.addParticle(mass=1)
    system.addParticle(mass=1)
    system.addParticle(mass=12)
    system.addParticle(mass=1)
    system.addParticle(mass=1)
    system.addParticle(mass=1)

    force = mm.NonbondedForce()
    force.setNonbondedMethod(force.NoCutoff)
    # These non-bonded parameters are completely wrong, but it does not matter
    # for the tests as long as we do not start testing the dynamic and
    # thermodynamics properties of methane.
    force.addParticle(charge=0, sigma=0.47, epsilon=3.5)
    force.addParticle(charge=0, sigma=0.47, epsilon=3.5)
    force.addParticle(charge=0, sigma=0.47, epsilon=3.5)
    force.addParticle(charge=0, sigma=0.47, epsilon=3.5)
    force.addParticle(charge=0, sigma=0.47, epsilon=3.5)
    force.addParticle(charge=0, sigma=0.47, epsilon=3.5)
    force.addParticle(charge=0, sigma=0.47, epsilon=3.5)
    force.addParticle(charge=0, sigma=0.47, epsilon=3.5)
    system.addForce(force)

    integrator = mm.LangevinIntegrator(300 * kelvin, 1 / picosecond, 2 * femtosecond)

    simulation = app.Simulation(topology, system, integrator)
    simulation.context.setPeriodicBoxVectors(*periodic_box_vector)
    simulation.context.setPositions(positions * nanometer)

    return simulation
