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
    Setup a minimal OpenMM simulation with two atoms.
    """
    periodic_box_vector = [
        [10,  0,  0],
        [ 0, 10,  0],
        [ 0,  0, 10],
    ]
    positions = np.array([[0, 0, 0], [0, 3, 0]], dtype=float)

    topology = app.Topology()
    chain = topology.addChain()
    element = app.Element.getBySymbol('C')
    residue = topology.addResidue(name='RES', chain=chain)
    topology.addAtom(element=element, name='A1', residue=residue)
    topology.addAtom(element=element, name='A2', residue=residue)

    system = mm.System()
    system.setDefaultPeriodicBoxVectors(*periodic_box_vector)
    system.addParticle(mass=72)
    system.addParticle(mass=72)

    force = mm.NonbondedForce()
    force.setNonbondedMethod(force.NoCutoff)
    force.addParticle(charge=0, sigma=0.47, epsilon=3.5)
    force.addParticle(charge=0, sigma=0.47, epsilon=3.5)
    system.addForce(force)

    integrator = mm.LangevinIntegrator(300 * kelvin, 1 / picosecond, 2 * femtosecond)

    simulation = app.Simulation(topology, system, integrator)
    simulation.context.setPeriodicBoxVectors(*periodic_box_vector)
    simulation.context.setPositions(positions * nanometer)

    return simulation
