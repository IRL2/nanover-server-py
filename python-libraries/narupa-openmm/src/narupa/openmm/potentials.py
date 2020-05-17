# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.
"""
Helpers around useful potentials that can
be applied to OpenMM systems.
"""

from typing import Iterable
from simtk import openmm as mm
from simtk.unit import Quantity, kilojoule_per_mole, nanometer


DEFAULT_RESTRAINT_FORCE_CONSTANT = 100 * kilojoule_per_mole / nanometer ** 2


def restraint_force(force_constant: Quantity = DEFAULT_RESTRAINT_FORCE_CONSTANT):
    """
    Generate an OpenMM force for position restraints.

    The position of the selected atoms is restrained with a harmonic potential
    the form of which is :math:`k * distance(r, r0)^2`; where :math:`k` is the
    force constant in kJ/(mol * nm^2), :math:`r` is the position of the atom,
    and :math:`r0` is the equilibrium position of that atom.

    >>> from simtk.unit import  kilojoule_per_mole, nanometer
    >>> # Create the force
    >>> force_constant = 225 * kilojoule_per_mole / nanometer ** 2
    >>> force = restraint_force(force_constant)
    >>> # Add particles to the force
    >>> force.addParticle(0, [10 * nanometer, 8 * nanometer, 2.3 * nanometer])
    >>> force.addParticle(2, [12 * nanometer, 8.4 * nanometer, 0.8 * nanometer])
    >>> # Add the force to the system
    >>> my_system.addForce(force)

    .. note::

        Like any force in OpenMM, the position restraints need to be added to
        the system before a context is created from the molecular system.

    :param force_constant: The force constant, :math:`k`. The value is expected
        to be a :class:`Quantity`, it will be converted to kJ/(mol * nm^2).
    :return: An OpenMM force to add to the system.
    """
    force_constant = force_constant.value_in_unit(kilojoule_per_mole / nanometer ** 2)
    force = mm.CustomExternalForce('k * periodicdistance(x, y, z, x0, y0, z0)^2')
    force.addGlobalParameter('k', force_constant)
    force.addPerParticleParameter('x0')
    force.addPerParticleParameter('y0')
    force.addPerParticleParameter('z0')
    return force


def restrain_particles(
        positions: Quantity,
        particle_indices: Iterable[int],
        force_constant: Quantity = DEFAULT_RESTRAINT_FORCE_CONSTANT,
):
    """
    Generate an OpenMM force that restrains the position of the selected particles.

    Apply a harmonic potential to restrain the position of the selected
    particles. The particles are restrained to their positions in the given
    position array (usually the initial positions).

    See :fun:`restraint_force` for the details of the potential.

    .. seealso: restraint_force

    :param positions: The positions of all the particles in the system.
    :param particle_indices: The indices in the system of the particles to
        restrain.
    :param force_constant: The force constant for the harmonic potential.
    :return: The populated OpenMM force to add to the OpenMM system.
    """
    force = restraint_force(force_constant)
    for index in particle_indices:
        x0, y0, z0 = positions[index].value_in_unit(nanometer)
        force.addParticle(index, [x0, y0, z0])
    return force

