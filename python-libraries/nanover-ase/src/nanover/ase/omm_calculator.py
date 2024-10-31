"""
ASE calculator for use with OpenMM.
"""

from typing import Optional

import nanover.openmm.serializer as serializer
import numpy as np
from ase import Atoms, Atom  # type: ignore
from ase.calculators.calculator import Calculator, all_changes
from openmm import System, State
from openmm.app import Topology, Simulation
from openmm.unit import angstrom, kilojoules_per_mole, kilojoule_per_mole, amu, Quantity

from nanover.ase.converter import KJMOL_TO_EV, ase_to_frame_data
from nanover.openmm.converter import add_openmm_topology_to_frame_data
from nanover.trajectory import FrameData


class OpenMMCalculator(Calculator):
    """
    Simple implementation of a ASE calculator for OpenMM. Initialises an OpenMM
    context with the given OpenMM simulation.

    :param simulation: An OpenMM simulation.
    :param atoms: ASE :class:`Atoms` to use with the calculator. The topology
        of the ASE atoms should be consistent with the OpenMM simulation.
        See :func:`~OpenMMCalculator.generate_atoms` for a helper function to
        generate a compatible ASE atoms object.

    """

    simulation: Simulation
    implemented_properties = ["energy", "forces"]

    def __init__(self, simulation, atoms: Optional[Atoms] = None, **kwargs):
        super().__init__(**kwargs)
        self.simulation = simulation
        self.context = self.simulation.context
        self.atoms = atoms

    @classmethod
    def from_xml(cls, input_xml, atoms: Optional[Atoms] = None, **kwargs):
        """
        Initialises an :class: OpenMMCalculator from a simulation serialised to
        XML with :module serializer.

        :param input_xml: XML file from which to create OpenMM simulation.
        :param atoms: ASE :class:`Atoms` to pass to the resulting OpenMMCalculator.
        :param kwargs: Keyword arguments for the OpenMMCalculator to be passed
            upon construction.
        :return: An :class: OpenMMCalculator.
        """
        with open(input_xml) as infile:
            simulation = serializer.deserialize_simulation(infile.read())
        return OpenMMCalculator(simulation, atoms, **kwargs)

    def calculate(
        self,
        atoms: Optional[Atoms] = None,
        properties=("energy", "forces"),
        system_changes=all_changes,
    ):
        if atoms is None:
            atoms = self.atoms
        if atoms is None:
            raise ValueError(
                "No ASE atoms supplied to calculator, and no ASE atoms supplied with initialisation."
            )

        self._set_positions(atoms.positions)
        energy, forces = self._calculate_openmm()
        self.results["energy"] = energy
        self.results["forces"] = forces

    def _calculate_openmm(self):
        state: State = self.context.getState(getEnergy=True, getForces=True)
        energy_kj_mol = state.getPotentialEnergy()
        energy = energy_kj_mol.value_in_unit(kilojoules_per_mole) * KJMOL_TO_EV
        forces_openmm = state.getForces(asNumpy=True)
        forces_angstrom = forces_openmm.value_in_unit(kilojoule_per_mole / angstrom)
        forces = forces_angstrom * KJMOL_TO_EV
        return energy, forces

    def _set_positions(self, positions):
        self.context.setPositions(positions * angstrom)

    def generate_atoms(self) -> Atoms:
        """
        Generates ASE atoms representation of the OpenMM system.

        :return: ASE :class:`Atoms`, with positions and chemical symbols set as
            according to the current state of the OpenMM system.
        """
        top: Topology = self.simulation.topology
        atoms = Atoms()
        system: System = self.simulation.system
        self.set_periodic_bounds(atoms, system)
        positions_openmm = self.context.getState(getPositions=True).getPositions(
            asNumpy=True
        )
        positions = positions_openmm.value_in_unit(angstrom)
        for openmm_atom in top.atoms():
            index = openmm_atom.index
            pos = positions[index]
            ase_atom = Atom(
                symbol=openmm_atom.element.symbol,
                position=pos,
                mass=openmm_atom.element.mass.value_in_unit(amu),
            )
            atoms.append(ase_atom)

        return atoms

    @property
    def topology(self):
        return self.simulation.topology

    @staticmethod
    def set_periodic_bounds(atoms: Atoms, system: System):
        """
        Sets ASE atoms object with the same periodic boundaries as that used in
        the given OpenMM system.

        :param atoms: ASE Atoms
        :param system: OpenMM system.
        """
        boxvectors: Quantity = system.getDefaultPeriodicBoxVectors()
        atoms.set_pbc(system.usesPeriodicBoundaryConditions())
        atoms.set_cell(
            np.array([vector.value_in_unit(angstrom) for vector in boxvectors])
        )


def openmm_ase_atoms_to_frame_data(
    ase_atoms: Atoms,
    *,
    topology: bool,
    **kwargs,
) -> FrameData:
    frame_data = ase_to_frame_data(
        ase_atoms,
        topology=False,
        **kwargs,
    )

    if topology:
        imd_calculator = ase_atoms.calc
        topology = imd_calculator.calculator.topology
        add_openmm_topology_to_frame_data(frame_data, topology)

    return frame_data
