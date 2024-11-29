"""
Calculate thermodynamic quantities associated with the simulation.
"""

from ase import Atoms, units


def compute_instantaneous_temperature(kinetic_energy: float, dof: int) -> float:
    """
    Calculate the instantaneous temperature of the system, using the same procedure as
    the :meth:`get_temperature` method of the :class:`Atoms` class in ASE. This
    function calculates the instantaneous temperature of the system from the kinetic
    energy of the system using the equation

    .. math::
        T = \frac{2 * \mathrm{KE}}{N_{\mathrm{dof}} * k_{\mathrm{B}}}

    where KE is the kinetic energy of the system (in eV), :math:`N_{dof}` is the number
    of degrees of freedom of the system and :math:`k_{\mathrm{B}}` is the Boltzmann
    constant in native ASE units.
    """

    return (2 * kinetic_energy) / (dof * units.kB)


def compute_dof(atoms: Atoms) -> int:
    """
    Compute the number of degrees of freedom of the system, using the same procedure as
    the :meth:`get_temperature` method of the :class:`Atoms` class in ASE.
    """
    dof = len(atoms) * 3
    for constraint in atoms._constraints:
        dof -= constraint.get_removed_dof(atoms)

    return dof
