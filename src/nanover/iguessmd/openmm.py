import os.path
import warnings
import xml.etree.ElementTree as ET
from os import PathLike
from typing import Any

import numpy as np
import numpy.typing as npt
from openmm import CustomCentroidBondForce, CustomExternalForce, OpenMMException
from openmm.app import Simulation

from nanover.iguessmd.base import iGUESSMDSimulation
from nanover.iguessmd.utils import (
    calculate_com,
    calculate_unit_tangents,
    calculate_weighted_centre,
    get_every_nth,
)
from nanover.openmm import serializer

iGUESSMD_FORCE_CONSTANT_PARAMETER_NAME = "iguessmd_k"
iGUESSMD_FORCE_CONSTANT_PARALLEL_PARAMETER_NAME = (
    f"{iGUESSMD_FORCE_CONSTANT_PARAMETER_NAME}_par"
)
iGUESSMD_FORCE_CONSTANT_PERPENDICULAR_PARAMETER_NAME = (
    f"{iGUESSMD_FORCE_CONSTANT_PARAMETER_NAME}_perp"
)
iGUESSMD_FORCE_EXPRESSION_SPHERICAL = (
    f"0.5 * {iGUESSMD_FORCE_CONSTANT_PARAMETER_NAME} * (dx^2 + dy^2 + dz^2)"
)
# tx, ty, tz are components of the unit tangent to the iGUESSMD path at the point x0, y0, z0
iGUESSMD_FORCE_EXPRESSION_PARALLEL = f"0.5 * {iGUESSMD_FORCE_CONSTANT_PARALLEL_PARAMETER_NAME} * (tx*dx + ty*dy + tz*dz)^2"
iGUESSMD_FORCE_EXPRESSION_PERPENDICULAR = f"0.5 * {iGUESSMD_FORCE_CONSTANT_PERPENDICULAR_PARAMETER_NAME} * ((dx^2 + dy^2 + dz^2) - (tx*dx + ty*dy + tz*dz)^2)"
iGUESSMD_FORCE_EXPRESSION_ATOM_NONPERIODIC = f"{iGUESSMD_FORCE_EXPRESSION_PARALLEL} + {iGUESSMD_FORCE_EXPRESSION_PERPENDICULAR}; dx=(x-x0); dy=(y-y0); dz=(z-z0)"
iGUESSMD_FORCE_EXPRESSION_ATOM_PERIODIC = f"{iGUESSMD_FORCE_EXPRESSION_PARALLEL} + {iGUESSMD_FORCE_EXPRESSION_PERPENDICULAR}; dx=(raw_dx - Lx*floor((raw_dx/Lx) + 0.5)); raw_dx=(x-x0); dy=(raw_dy - Ly*floor((raw_dy/Ly) + 0.5)); raw_dy=(y-y0); dz=(raw_dz - Lz*floor((raw_dz/Lz) + 0.5)); raw_dz=(z-z0)"
iGUESSMD_FORCE_EXPRESSION_COM_NONPERIODIC = f"{iGUESSMD_FORCE_EXPRESSION_PARALLEL} + {iGUESSMD_FORCE_EXPRESSION_PERPENDICULAR}; dx=((x1-x2)-x0); dy=((y1-y2)-y0); dz=((z1-z2)-z0)"
iGUESSMD_FORCE_EXPRESSION_COM_PERIODIC = f"{iGUESSMD_FORCE_EXPRESSION_PARALLEL} + {iGUESSMD_FORCE_EXPRESSION_PERPENDICULAR}; dx=(raw_dx - Lx*floor((raw_dx/Lx) + 0.5)); raw_dx=((x1-x2)-x0); dy=(raw_dy - Ly*floor((raw_dy/Ly) + 0.5)); raw_dy=((y1-y2)-y0); dz=(raw_dz - Lz*floor((raw_dz/Lz) + 0.5)); raw_dz=((z1-z2)-z0)"


class OMMiGUESSMDSimulation(iGUESSMDSimulation):
    """
    A class for performing iGUESSMD on an OpenMM simulation.
    """

    @classmethod
    def from_simulation(
        cls,
        simulation: Simulation,
        iguessmd_mobile_atom_indices: np.ndarray,
        iguessmd_reaction_path: np.ndarray,
        iguessmd_force_constant: float | np.ndarray,
        *,
        iguessmd_ref_atom_indices: np.ndarray | None = None,
        iguessmd_coord_type: str | None = None,
        iguessmd_mobile_atom_weights: np.ndarray | None = None,
        iguessmd_ref_atom_weights: np.ndarray | None = None,
        name: str | None = None,
    ):
        """
        Construct the iGUESSMD simulation from an existing OpenMM simulation.

        :param simulation: An existing OpenMM Simulation
        :param iguessmd_mobile_atom_indices: A NumPy array of indices of the mobile atoms to which the
          iGUESSMD force should be applied (0-D or 1-D for single atom, 1-D for COM)
        :param iguessmd_reaction_path: A NumPy array of coordinates defining the path that the
          iGUESSMD force will take during the iGUESSMD simulation
        :param iguessmd_force_constant: The force constant(s) defining the iGUESSMD force
          (either a float for a spherically symmetric harmonic potential or a 1-D NumPy
           array defining a parallel and perpendicular force constant for an asymmetric
           harmonic potential)
        :param iguessmd_ref_atom_indices: A NumPy array of indices of the atoms defining the centre
           used to define the relative reaction path of the mobile atoms, to which the iGUESSMD
           is also applied if relative coordinates are used (0-D or 1-D for single atom,
           1-D for COM)
        :param iguessmd_coord_type: The type of coordinate system used to define the iGUESSMD
           reaction path.
           Possible values:
             - 'abs' : indicates that the reaction path is given in the same Cartesian coordinate
                       system of the simulation, with the origin as the reference point R_ref(r)
             - 'rel' : indicates that the reaction path is given in the relative coordinate
                       system defined by the centres R(r) and R_ref(r), which are defined
                       respectively by the mobile and reference atom positions
        :param iguessmd_mobile_atom_weights: [Optional] An array of weights used to define R(r)
           from the mobile atom positions (must sum to 1). If no weights are given, the weights
           are defined to yield R(r) as the centre of mass of the mobile atom group.
        :param iguessmd_ref_atom_weights: [Optional] An array of weights used to define R_ref(r)
           from the reference atom positions (must sum to 1). If no weights are given, the weights
           are defined to yield R(r) as the centre of mass of the reference atom group.
        :param name: An optional name for the simulation
        """

        # Instantiate class
        sim = cls(name)

        # Define OpenMM Simulation
        sim.simulation = simulation

        # Initialise all objects relevant to the iGUESSMD simulation
        sim._initialise_iguessmd_simulation(
            iguessmd_mobile_atom_indices,
            iguessmd_reaction_path,
            iguessmd_force_constant,
            iguessmd_ref_atom_indices,
            iguessmd_coord_type,
            iguessmd_mobile_atom_weights,
            iguessmd_ref_atom_weights,
        )

        return sim

    @classmethod
    def from_xml_path(
        cls,
        path: PathLike | str,
        iguessmd_mobile_atom_indices: np.ndarray,
        iguessmd_reaction_path: np.ndarray,
        iguessmd_force_constant: float | np.ndarray,
        *,
        iguessmd_ref_atom_indices: np.ndarray | None = None,
        iguessmd_coord_type: str | None = None,
        iguessmd_mobile_atom_weights: np.ndarray | None = None,
        iguessmd_ref_atom_weights: np.ndarray | None = None,
        name: str | None = None,
    ):
        """
        Construct the iGUESSMD simulation from an existing NanoVer OpenMM XML file located at a given path.

        :param path: Path to the NanoVer OpenMM XML file defining the simulation
        :param iguessmd_mobile_atom_indices: A NumPy array of indices of the atoms to which the
          iGUESSMD force should be applied (0-D or 1-D for single atom, 1-D for COM)
        :param iguessmd_reaction_path: A NumPy array of coordinates defining the path that the
          iGUESSMD force will take during the iGUESSMD simulation
        :param iguessmd_force_constant: The force constant(s) defining the iGUESSMD force
          (either a float for a spherically symmetric harmonic potential or a 1-D NumPy
           array defining a parallel and perpendicular force constant for an asymmetric
           harmonic potential)
        :param iguessmd_ref_atom_indices: A NumPy array of indices of the atoms defining the centre
           used to define the relative reaction path of the mobile atoms, to which the iGUESSMD
           is also applied if relative coordinates are used (0-D or 1-D for single atom,
           1-D for COM)
        :param iguessmd_coord_type: The type of coordinate system used to define the iGUESSMD
           reaction path.
           Possible values:
             - 'abs' : indicates that the reaction path is given in the same Cartesian coordinate
                       system of the simulation, with the origin as the reference point R_ref(r)
             - 'rel' : indicates that the reaction path is given in the relative coordinate
                       system defined by the centres R(r) and R_ref(r), which are defined
                       respectively by the mobile and reference atom positions
        :param iguessmd_mobile_atom_weights: [Optional] An array of weights used to define R(r)
           from the mobile atom positions (must sum to 1). If no weights are given, the weights
           are defined to yield R(r) as the centre of mass of the mobile atom group.
        :param iguessmd_ref_atom_weights: [Optional] An array of weights used to define R_ref(r)
           from the reference atom positions (must sum to 1). If no weights are given, the weights
           are defined to yield R(r) as the centre of mass of the reference atom group.
        :param name: An optional name for the simulation
        """

        # Instantiate class
        sim = cls(name)

        # Load OpenMM Simulation from the xml path
        sim.xml_path = path
        assert sim.xml_path is not None
        with open(sim.xml_path) as infile:
            sim.simulation = serializer.deserialize_simulation(infile)

        # Initialise all objects relevant to the iGUESSMD simulation
        sim._initialise_iguessmd_simulation(
            iguessmd_mobile_atom_indices,
            iguessmd_reaction_path,
            iguessmd_force_constant,
            iguessmd_ref_atom_indices,
            iguessmd_coord_type,
            iguessmd_mobile_atom_weights,
            iguessmd_ref_atom_weights,
        )

        return sim

    def __init__(self, name: str | None = None):
        # Inherit most iGUESSMD attributes from base class
        super().__init__(name)

        # Attributes from loading from file
        self.xml_path: PathLike | str | None = None

        # Attributes specific to OpenMM Simulation
        self.simulation: Simulation
        self.checkpoint: Any | None = None
        self.iguessmd_force: CustomCentroidBondForce | None = None
        self.iguessmd_force_groups: list
        self.loaded_iguessmd_force_from_sim: bool = False

        self._sim_uses_pbcs: bool | None = None
        self._pbc_box_lengths: np.ndarray | None = None

        self.n_iguessmd_atom_indices: int | None = None

    def _initialise_iguessmd_simulation(
        self,
        iguessmd_atom_indices: np.ndarray,
        iguessmd_path: np.ndarray,
        iguessmd_force_constant: float | np.ndarray,
        iguessmd_ref_atom_indices: np.ndarray | None = None,
        iguessmd_coord_type: str | None = None,
        iguessmd_mobile_atom_weights: np.ndarray | None = None,
        iguessmd_ref_atom_weights: np.ndarray | None = None,
    ):
        """
        Set the fields relevant to the iGUESSMD simulation and add the iGUESSMD force to the system.
        Called upon when constructing an OMMiGUESSMDSimulation using the classmethods
        from_simulation or from_xml_path.
        """
        assert iguessmd_atom_indices.size >= 1
        if iguessmd_atom_indices.size == 1 and iguessmd_atom_indices.shape == ():
            iguessmd_atom_indices = np.reshape(iguessmd_atom_indices, (1,))

        self.mobile_atom_indices = iguessmd_atom_indices
        if iguessmd_mobile_atom_weights is not None and isinstance(
            iguessmd_mobile_atom_weights, np.ndarray
        ):
            self._check_custom_weights_are_normalised(iguessmd_mobile_atom_weights)
            self.mobile_atom_weights = iguessmd_mobile_atom_weights

        # Define atom weights for COM weights if no weights passed
        if self.mobile_atom_indices is not None and self.mobile_atom_weights is None:
            self.mobile_atom_weights = self._calculate_com_weights(
                self.mobile_atom_indices
            )
        self.n_iguessmd_atom_indices = self.mobile_atom_indices.size
        self.reaction_path = iguessmd_path

        if iguessmd_ref_atom_indices is not None:
            assert iguessmd_ref_atom_indices.size >= 1
            if (
                iguessmd_ref_atom_indices.size == 1
                and iguessmd_ref_atom_indices.shape == ()
            ):
                iguessmd_ref_atom_indices = np.reshape(iguessmd_ref_atom_indices, (1,))
            self.ref_atom_indices = iguessmd_ref_atom_indices
            if iguessmd_ref_atom_weights is not None:
                self._check_custom_weights_are_normalised(iguessmd_ref_atom_weights)
                self.ref_atom_weights = iguessmd_ref_atom_weights

            if self.ref_atom_indices is not None and self.ref_atom_weights is None:
                self.ref_atom_weights = self._calculate_com_weights(
                    self.ref_atom_indices
                )

        # Allow user to choose iGUESSMD coordinate type (absolute or relative) with input
        if iguessmd_coord_type is not None:
            self._coord_type = iguessmd_coord_type
        # If no iGUESSMD coordinate type specified, use relative coordinates if
        # reference group is specified
        elif self.ref_atom_indices is not None:
            self._coord_type = "rel"
        # Otherwise, default to absolute coordinates
        else:
            self._coord_type = "abs"

        # Check if simulation employs periodic boundary conditions
        self._sim_uses_pbcs = self.simulation.system.usesPeriodicBoundaryConditions()
        if self._sim_uses_pbcs:
            # TODO: only works with orthorhombic PBCs at the moment, need to generalise
            self._pbc_box_lengths = np.diag(
                self.simulation.context.getState()
                .getPeriodicBoxVectors(asNumpy=True)
                ._value
            )

        if (
            isinstance(iguessmd_force_constant, np.ndarray)
            and iguessmd_force_constant.size == 2
        ):
            self.k = iguessmd_force_constant
            self.k_par = iguessmd_force_constant[0]
            self.k_perp = iguessmd_force_constant[1]
        else:
            self.k = self.k_par = self.k_perp = iguessmd_force_constant

        self.define_iguessmd_simulation_atom_positions_arrays()

        self.current_iguessmd_force_position = self.reaction_path[0]
        self.current_iguessmd_force_position_index = 0

        self.calculate_iguessmd_path_tangents()
        self.current_iguessmd_force_tangent = self.reaction_path_unit_tangents[0]

        # Check whether iGUESSMD force is already present
        self.check_for_existing_iguessmd_force()

        if not self.loaded_iguessmd_force_from_sim:
            # Create iGUESSMD force and add it to the system
            self.add_iguessmd_force_to_system()

        assert self.iguessmd_force is not None

        self.simulation.context.reinitialize(preserveState=True)

        # Create a checkpoint of the simulation
        self.checkpoint = self.simulation.context.createCheckpoint()

    def reset(self):
        """
        Reset the iGUESSMD simulation to its initial state, and reset the arrays output
        by the iGUESSMD simulation.
        """
        assert (
            self.simulation is not None
            and self.checkpoint is not None
            and self.reaction_path is not None
            and self.iguessmd_force is not None
        )

        # Reset iGUESSMD force position
        self.current_iguessmd_force_position = self.reaction_path[0]
        self.current_iguessmd_force_position_index = 0
        self.update_iguessmd_force_position()

        # Reset simulation, reinitialise context to be safe
        self.simulation.context.reinitialize()
        self.simulation.context.loadCheckpoint(self.checkpoint)

        # Reset or remove iGUESSMD simulation arrays
        self.define_iguessmd_simulation_atom_positions_arrays()
        del self.iguessmd_sim_forces
        del self.iguessmd_sim_work_done

    def remove_iguessmd_force_from_system(self):
        """
        Remove any iGUESSMD forces from the system.
        """
        forces = self.simulation.system.getForces()
        forces_to_remove = []
        for i in range(len(forces)):
            if type(forces[i]) is type(self.iguessmd_force):
                try:
                    if (
                        forces[i].getGlobalParameterName(0)
                        == iGUESSMD_FORCE_CONSTANT_PARALLEL_PARAMETER_NAME
                    ):
                        forces_to_remove.append(i)
                except OpenMMException:
                    continue

        # Remove any iGUESSMD forces, accounting for the changes in indices
        # as forces are removed
        for j in range(len(forces_to_remove)):
            self.simulation.system.removeForce(forces_to_remove[j] - j)

    def calculate_iguessmd_path_tangents(self):
        """
        Calculate the unit tangent vectors of the iGUESSMD path using the forward
        difference approximation.
        """

        # Check that path exists
        assert self.reaction_path is not None

        # Calculate unit tangents
        self.reaction_path_unit_tangents = calculate_unit_tangents(self.reaction_path)

    def get_iguessmd_atom_positions(self):
        """
        Retrieve the positions of the atoms with which the iGUESSMD force is
        interacting, and add them to the array of positions to save.
        """
        positions = self.simulation.context.getState(
            getPositions=True, enforcePeriodicBox=False
        ).getPositions(asNumpy=True)
        self.iguessmd_sim_mobile_atom_positions[
            self.current_iguessmd_force_position_index
        ] = positions[self.mobile_atom_indices]
        if (
            self._coord_type == "rel"
            and self.ref_atom_indices is not None
            and self.iguessmd_sim_ref_atom_positions is not None
        ):
            self.iguessmd_sim_ref_atom_positions[
                self.current_iguessmd_force_position_index
            ] = positions[self.ref_atom_indices]

    def run_equilibration_with_initial_restraint(self, n_steps: int):
        """
        Perform an equilibration of the system with the restraint fixed at the
        initial position defined in the iGUESSMD force path.

        :param n_steps: Number of simulation steps to run equilibration for.
        """

        try:
            assert np.all(self.current_iguessmd_force_position == self.reaction_path[0])
        except AssertionError:
            raise AssertionError(
                "Restraint is not located at the initial position "
                "of the iGUESSMD force, equilibration aborted. To prepare the system "
                "appropriately, the restraint should be placed at the first point "
                "on the path that the iGUESSMD force will take during the iGUESSMD simulation."
            )

        self.simulation.step(n_steps)

    def save_simulation(
        self,
        output_filepath: PathLike | str,
        save_state: bool = False,
        save_iguessmd_force: bool = False,
    ):
        """
        Save the simulation to a NanoVer OpenMM XML file, with the option to include the
        iGUESSMD force in the XML file.

        :param output_filepath: Path to output file to save the simulation to.
        :param save_state: If True, save the present state of the simulation to the XML file (Optional).
        :param save_iguessmd_force: Bool defining whether to save the iGUESSMD force in the XML file (Optional).
        """
        assert output_filepath is not None

        if save_iguessmd_force:
            with open(output_filepath, "w") as outfile:
                outfile.write(
                    serializer.serialize_simulation(
                        self.simulation, save_state=save_state
                    )
                )

        else:
            # Temporarily remove iGUESSMD force from simulation
            self.remove_iguessmd_force_from_system()

            # Save simulation without SMD force
            with open(output_filepath, "w") as outfile:
                xml_string = serializer.serialize_simulation(
                    self.simulation, save_state=save_state
                )

                # Manually remove parameters for now...
                params_to_remove = {
                    iGUESSMD_FORCE_CONSTANT_PARAMETER_NAME,
                    iGUESSMD_FORCE_CONSTANT_PARALLEL_PARAMETER_NAME,
                    iGUESSMD_FORCE_CONSTANT_PERPENDICULAR_PARAMETER_NAME,
                    *({"Lx", "Ly", "Lz"} if self._sim_uses_pbcs else set()),
                    *({"x2", "y2", "z2"} if self._coord_type == "abs" else set()),
                }

                xml_root = ET.fromstring(xml_string)

                for elem in xml_root.iter():
                    for param in params_to_remove:
                        elem.attrib.pop(param, None)

                # Regenerate xml string and write to file
                xml_string = ET.tostring(xml_root, encoding="unicode")
                outfile.write(xml_string)

            # Add iGUESSMD force back to the system
            self.add_iguessmd_force_to_system()

    def generate_starting_structures(
        self,
        interval_ps: float,
        n_structures: int,
        output_directory: PathLike | str | None = None,
        filename_prefix: str | None = None,
        save_iguessmd_force: bool = False,
    ):
        """
        Generate the specified number of starting structures by running the simulation for the specified
        interval with the initial restraint applied to the system and saving structures at regular intervals.
        Structures are saved to the output path, optionally with the filename prefix. If no output path is
        specified, structures will be saved to the current working directory. If no prefix is specified,
        the output files will be named automatically.

        :param interval_ps: Interval to run the simulation for, in picoseconds.
        :param n_structures: Number of structures to generate.
        :param output_directory: Output directory to save the structures to (Optional).
        :param filename_prefix: Prefix for output files (Optional).
        :param save_iguessmd_force: Bool defining whether to save the iGUESSMD force in the XML file (Optional).
        """
        timestep_ps = self.simulation.integrator.getStepSize()._value
        n_steps_struct_interval = int(
            np.floor(interval_ps / (n_structures * timestep_ps))
        )

        if output_directory is None:
            output_directory = os.getcwd()
            print(
                f"No output directory specified. Generated starting structures "
                f"will be saved to {output_directory}."
            )

        if filename_prefix is None:
            filename_prefix = "iguessmd_structure"
            print(
                f"No filename prefix specified. Generated starting structures "
                f"will be saved with the prefix {filename_prefix}."
            )

        print(
            f"Generating {n_structures} structures in {interval_ps} ps simulation...\n"
            f"Structures will be saved to {output_directory}\n"
        )

        for i in range(n_structures):
            # Run set of simulation steps to generate next structure
            self.simulation.step(n_steps_struct_interval)

            # Save structure to output file
            outfile_path = os.path.join(
                output_directory, filename_prefix + "_" + str(i + 1) + ".xml"
            )
            self.save_simulation(
                output_filepath=outfile_path,
                save_state=True,
                save_iguessmd_force=save_iguessmd_force,
            )

        print(f"Structure generation complete: {n_structures} structures generated.")

    def run_iguessmd(self, progress_interval: int = 100000):
        """
        Perform an iGUESSMD simulation on the system using the iGUESSMD force and
        path defined, store the positions of the atoms with which the
        iGUESSMD force interacts, and calculate the work done along the reaction
        coordinate.

        :param progress_interval: Interval defining how regularly to print the progress
          of the simulation, in terms of the number of simulation steps.
        """

        # Get initial atom positions for iGUESSMD force at initial position
        assert (
            self.iguessmd_sim_mobile_atom_positions is not None
            and self.iguessmd_force is not None
            and np.all(self.current_iguessmd_force_position == self.reaction_path[0])
            and self.current_iguessmd_force_position_index == 0
        )
        self.get_iguessmd_atom_positions()

        # Run iGUESSMD procedure
        n_steps = self.reaction_path.shape[0]

        print("Starting SMD simulation...")
        for step in range(1, n_steps):
            # Perform single simulation step
            self.simulation.step(1)

            # Update force position index
            self.current_iguessmd_force_position_index = step

            # Update iGUESSMD force position
            self.update_iguessmd_force_position()

            # Retrieve atom positions on step
            self.get_iguessmd_atom_positions()

            # Print step at intervals
            if step % progress_interval == 0:
                print(f"Steps completed: {step}")

        print("iGUESSMD simulation completed. Calculating work done...")

        # Calculate the work done along the reaction coordinate
        self.calculate_cumulative_work_done()

        print("Work done calculated.")

    def _calculate_iguessmd_forces(self, interaction_centre_positions):
        """
        Calculate the iGUESSMD forces that acted on the system during the simulation in kJ mol-1 nm-1.

        :param interaction_centre_positions: Array of positions defining the centre
          (single atom or COM of group of atoms) with which the iGUESSMD force interacted during the simulation.
        """

        assert (
            interaction_centre_positions is not None
            and np.all(self.reaction_path.shape == interaction_centre_positions.shape)
            and self.k_par is not None
            and self.k_perp is not None
        )

        # Calculate displacements and tangential components
        displacements = interaction_centre_positions - self.reaction_path
        tangential_displacements = self._calculate_tangential_displacements(
            displacements
        )

        # Calculate force component along RC
        parallel_forces = -self.k_par * tangential_displacements
        perpendicular_forces = -self.k_perp * (displacements - tangential_displacements)

        self.iguessmd_sim_forces = parallel_forces + perpendicular_forces

    def _calculate_tangential_displacements(
        self, displacements: np.ndarray
    ) -> np.ndarray:
        """
        Calculate the tangential components of the displacement vectors between
        the iGUESSMD force origin and the centre to which the iGUESSMD force is
        attached, given an array of displacements.
        :param displacements: NumPy array of displacement vectors
        :return: NumPy array of tangential displacements
        """
        assert self.reaction_path_unit_tangents is not None

        tangential_displacements = (
            (
                np.linalg.vecdot(
                    displacements, self.reaction_path_unit_tangents, axis=1
                )[:, None]
            )
            * self.reaction_path_unit_tangents
        )

        return tangential_displacements

    def _calculate_work_done(self):
        """
        Calculate the work done along the reaction coordinate in kJ mol-1.
        """
        assert self.reaction_path is not None and self.iguessmd_sim_forces is not None
        smd_force_displacements = np.diff(self.reaction_path, axis=0)
        work_done_array = np.zeros(self.iguessmd_sim_forces.shape[0])
        work_done_array[1:] = np.linalg.vecdot(
            self.iguessmd_sim_forces[:-1], smd_force_displacements, axis=1
        )
        self.iguessmd_sim_work_done = np.cumsum(work_done_array, axis=0)

    def save_iguessmd_simulation_data(
        self,
        path: PathLike | str | None = None,
        every_nth: int | None = None,
        include_end: bool = False,
        save_work_done: bool = True,
        save_x_rel_positions: bool = True,
        save_centre_positions: bool = True,
        save_atom_positions: bool = True,
        work_done_dtype: npt.DTypeLike = np.float32,
        x_rel_positions_dtype: npt.DTypeLike = np.float32,
        centre_positions_dtype: npt.DTypeLike = np.float32,
        atom_positions_dtype: npt.DTypeLike = np.float32,
    ):
        """
        Save the data produced by the iGUESSMD simulation in binary form that can be read
        into NumPy arrays. The following data are saved, in the order listed below:

        - Trajectories of the atoms to which the iGUESSMD force was applied, in nm
        - Work done along the reaction coordinate defined by the path of the iGUESSMD force, in kJ mol-1

        :param path: Path to the file to which the data will be saved
        :param every_nth:  (int | None) only save the values of work and/or positions at every
          nth point along the trajectory
        :param include_end: (Bool) whether to save final values of the work and/or positions
          regardless of stride defined by every nth
        :param save_work_done: Bool determining whether to save the work done
        :param save_x_rel_positions: Bool determining whether to save the positions of the x_rel
          coordinate defined during the iGUESSMD simulation
        :param save_centre_positions: Bool determining whether to save the positions of the centre R(r)
          defined by the mobile atoms (and the centre R_ref(r) defined by the reference atoms if relative
          coordinates used)
        :param save_atom_positions: Bool determining whether to save the positions of the atom(s)
        :param work_done_dtype: Data type of the work done array to save
        :param x_rel_positions_dtype: Data type of the x_rel positions array to save
        :param centre_positions_dtype: Data type of the centre positions array(s) to save
        :param atom_positions_dtype: Data type of the atom positions array(s) to save
        """

        if path is None:
            raise ValueError(
                "Output file path cannot be None. Please specify an output file path."
            )

        if self.iguessmd_sim_work_done is None:
            raise ValueError(
                "Missing values for the work done. This data can only be saved after "
                "the iGUESSMD calculation is completed."
            )

        if self.iguessmd_sim_x_rel_positions is None:
            raise ValueError(
                "Missing values for the x_rel positions. This data can only be saved after "
                "the iGUESSMD calculation is completed."
            )

        if self.iguessmd_sim_mobile_atom_positions is None or np.all(
            self.iguessmd_sim_mobile_atom_positions == 0.0
        ):
            raise ValueError(
                "Missing values for the mobile atom positions. This data can only be saved after "
                "the iGUESSMD calculation is completed."
            )

        if self.iguessmd_sim_mobile_centre_positions is None or np.all(
            self.iguessmd_sim_mobile_centre_positions == 0.0
        ):
            raise ValueError(
                "Missing values for the mobile centre positions. This data can only be saved after"
                "the iGUESSMD calculation is completed."
            )

        if self._coord_type == "rel":
            if self.iguessmd_sim_ref_atom_positions is None or np.all(
                self.iguessmd_sim_ref_atom_positions == 0.0
            ):
                raise ValueError(
                    "Missing values for the reference atom positions. This data can only be saved after "
                    "the iGUESSMD calculation is completed."
                )

            if self.iguessmd_sim_ref_centre_positions is None or np.all(
                self.iguessmd_sim_ref_centre_positions == 0.0
            ):
                raise ValueError(
                    "Missing values for the reference centre positions. This data can only be saved after"
                    "the iGUESSMD calculation is completed."
                )

        if every_nth is None:
            every_nth = 1

        if (
            every_nth is not None
            and include_end
            and (self.iguessmd_sim_work_done.size - 1) % every_nth != 0
        ):
            warnings.warn(
                "Choice of every_nth yields different time step between the "
                "final two array entries compared to the rest of the trajectory."
            )

        iguessmd_simulation_data: dict[str, Any] = {}
        iguessmd_simulation_data["data_timestep_ps"] = (
            every_nth * self.simulation.integrator.getStepSize()._value
        )

        # Optionally save work done
        if save_work_done:
            iguessmd_simulation_data["iguessmd_sim_work_done"] = get_every_nth(
                self.iguessmd_sim_work_done.astype(work_done_dtype),
                axis=0,
                every_nth=every_nth,
                include_end=include_end,
            )
            print("Work done saved to simulation data file.")

        # Optionally save the x_rel positions
        if save_x_rel_positions:
            iguessmd_simulation_data["iguessmd_sim_x_rel_positions"] = get_every_nth(
                self.iguessmd_sim_x_rel_positions,
                axis=0,
                every_nth=every_nth,
                include_end=include_end,
            ).astype(x_rel_positions_dtype)

        # Optionally save atom positions
        if save_atom_positions:
            iguessmd_simulation_data["iguessmd_sim_mobile_atom_positions"] = (
                get_every_nth(
                    self.iguessmd_sim_mobile_atom_positions,
                    axis=0,
                    every_nth=every_nth,
                    include_end=include_end,
                ).astype(atom_positions_dtype)
            )
            if self._coord_type == "rel":
                assert self.iguessmd_sim_ref_atom_positions is not None
                iguessmd_simulation_data["iguessmd_sim_ref_atom_positions"] = (
                    get_every_nth(
                        self.iguessmd_sim_ref_atom_positions,
                        axis=0,
                        every_nth=every_nth,
                        include_end=include_end,
                    ).astype(atom_positions_dtype)
                )
            print("Atom positions saved to simulation data file.")

        if save_centre_positions:
            iguessmd_simulation_data["iguessmd_sim_mobile_centre_positions"] = (
                get_every_nth(
                    self.iguessmd_sim_mobile_centre_positions.astype(
                        centre_positions_dtype
                    ),
                    axis=0,
                    every_nth=every_nth,
                    include_end=include_end,
                )
            )
            if self._coord_type == "rel":
                assert self.iguessmd_sim_ref_centre_positions is not None
                iguessmd_simulation_data["iguessmd_sim_ref_centre_positions"] = (
                    get_every_nth(
                        self.iguessmd_sim_ref_centre_positions.astype(
                            centre_positions_dtype
                        ),
                        axis=0,
                        every_nth=every_nth,
                        include_end=include_end,
                    )
                )
            print("Centre positions saved to simulation data file.")

        np.savez_compressed(path, **iguessmd_simulation_data)

    def save_general_iguessmd_data(self, path: PathLike | str | None = None):
        """
        Saves general data related to the iGUESSMD simulation in binary form that can be read
        into NumPy arrays. The following data are saved, in the order listed below:

        - Indices of the atoms to which the iGUESSMD force is applied
        - Positions defining the path of the iGUESSMD force, in nm
        - The type of coordinates used to define the iGUESSMD reaction path
        - Force constant of the iGUESSMD force, in kJ mol-1 nm-2
        - Temperature of the simulation, in Kelvin
        - Time step of the simulation, in picoseconds

        :param path: Path to the file to which the data will be saved.
        """

        if path is None:
            raise ValueError(
                "Output file path cannot be None. Please specify an output file path."
            )

        # Create dictionary containing general iGUESSMD data
        general_iguessmd_data: dict[str, Any] = {
            "mobile_atom_indices": self.mobile_atom_indices,
            "reaction_path": self.reaction_path,
            "x_rel_coordinate_type": self._coord_type,
            "force_constant_kJ/(mol*nm2)": self.k,
            "temperature_K": self.simulation.integrator.getTemperature()._value,
            "timestep_ps": self.simulation.integrator.getStepSize()._value,
        }

        # Add the indices of the reference atoms (if relative coordinate used)
        if self._coord_type == "rel" and self.ref_atom_indices is not None:
            general_iguessmd_data["ref_atom_indices"] = self.ref_atom_indices

        # Save dictionary to file defined by path
        np.savez_compressed(path, **general_iguessmd_data)

    def check_for_existing_iguessmd_force(self):
        """
        Check whether the loaded simulation already contains an iGUESSMD force of
        the correct type.
        """
        try:
            force_constant_parallel = self.simulation.context.getParameter(
                iGUESSMD_FORCE_CONSTANT_PARALLEL_PARAMETER_NAME
            )
            force_constant_perpendicular = self.simulation.context.getParameter(
                iGUESSMD_FORCE_CONSTANT_PERPENDICULAR_PARAMETER_NAME
            )
            n_forces = self.simulation.system.getNumForces()
            iguessmd_force = self.simulation.system.getForce(n_forces - 1)
            assert isinstance(iguessmd_force, CustomCentroidBondForce)
            assert (
                iguessmd_force.getEnergyFunction()
                == iGUESSMD_FORCE_EXPRESSION_COM_NONPERIODIC
                or iguessmd_force.getEnergyFunction()
                == iGUESSMD_FORCE_EXPRESSION_COM_PERIODIC
            )
            # TODO: Think about how to deal with multiple centres cleanly
            n_force_groups = iguessmd_force.getNumGroups()
            assert n_force_groups == 2 if self.ref_atom_indices is not None else 1
            assert force_constant_parallel == self.k_par
            assert force_constant_perpendicular == self.k_perp
            assert np.array_equal(
                np.asarray(iguessmd_force.getGroupParameters(0)[0]), self.mobile_atom_indices
            )
            assert np.array_equal(
                np.asarray(iguessmd_force.getGroupParameters(0)[1]), self.mobile_atom_weights
            )
            if n_force_groups == 2:
                assert np.array_equal(
                    np.asarray(iguessmd_force.getGroupParameters(1)[0]), self.ref_atom_indices
                )
                assert np.array_equal(
                    np.asarray(iguessmd_force.getGroupParameters(1)[1]), self.ref_atom_weights
                )
            assert np.all(
                iguessmd_force.getBondParameters(0)[1]
                == np.concatenate(
                    [
                        self.current_iguessmd_force_position,
                        self.current_iguessmd_force_tangent,
                    ]
                )
            )
            print("iGUESSMD force already present in loaded simulation.")
            self.iguessmd_force = iguessmd_force
            self.loaded_iguessmd_force_from_sim = True
            self.iguessmd_force_groups = [i for i in range(n_force_groups)]

        except OpenMMException:
            self.loaded_iguessmd_force_from_sim = False

    def add_iguessmd_force_to_system(self):
        assert (
            self.current_iguessmd_force_position is not None
            and self.current_iguessmd_force_tangent is not None
            and self.k_par is not None
            and self.k_perp is not None
            and self._sim_uses_pbcs is not None
            and self._coord_type is not None
        )

        x0, y0, z0 = self.current_iguessmd_force_position
        tx, ty, tz = self.current_iguessmd_force_tangent
        iguessmd_force = iguessmd_com_force(
            self.k_par,
            self.k_perp,
            self._sim_uses_pbcs,
            self._pbc_box_lengths,
            self._coord_type,
        )
        force_groups = [
            iguessmd_force.addGroup(self.mobile_atom_indices, self.mobile_atom_weights)
        ]
        if self._coord_type == "rel":
            force_groups.append(
                iguessmd_force.addGroup(self.ref_atom_indices, self.ref_atom_weights)
            )
        self.iguessmd_force_groups = force_groups
        iguessmd_force.addBond(force_groups, [x0, y0, z0, tx, ty, tz])
        self.iguessmd_force = iguessmd_force
        self.simulation.system.addForce(self.iguessmd_force)
        self.simulation.context.reinitialize(preserveState=True)

    def update_iguessmd_force_position(self):
        self.current_iguessmd_force_position = self.reaction_path[
            self.current_iguessmd_force_position_index
        ]
        self.current_iguessmd_force_tangent = self.reaction_path_unit_tangents[
            self.current_iguessmd_force_position_index
        ]
        assert (
            self.current_iguessmd_force_position is not None
            and self.current_iguessmd_force_tangent is not None
            and self.iguessmd_force is not None
            and self.iguessmd_force_groups is not None
        )
        x0, y0, z0 = self.current_iguessmd_force_position
        tx, ty, tz = self.current_iguessmd_force_tangent
        self.iguessmd_force.setBondParameters(
            0, self.iguessmd_force_groups, [x0, y0, z0, tx, ty, tz]
        )
        self.iguessmd_force.updateParametersInContext(self.simulation.context)

    def define_iguessmd_simulation_atom_positions_arrays(self):
        assert self.mobile_atom_indices is not None and self.reaction_path is not None
        self.iguessmd_sim_mobile_atom_positions = np.zeros(
            (self.reaction_path.shape[0], self.mobile_atom_indices.size, 3)
        )
        if self._coord_type == "rel" and self.ref_atom_indices is not None:
            self.iguessmd_sim_ref_atom_positions = np.zeros(
                (self.reaction_path.shape[0], self.ref_atom_indices.size, 3)
            )

    def _calculate_com_weights(self, atom_indices: np.ndarray) -> np.ndarray:
        """
        Calculate the relative weights of the atomic positions defining the centre of mass of the
        mobile (and reference if available) atom groups.
        :param atom_indices: NumPy array containing the indices of the atoms defining the group
          for which to calculate the mass weights.
        :return: NumPy array containing the mass weights of the atoms in the selection
        """
        assert atom_indices is not None and self.simulation is not None
        atom_masses = np.zeros(atom_indices.size)
        for index in range(atom_indices.size):
            atom_masses[index] = self.simulation.system.getParticleMass(
                atom_indices[index]
            )._value
        atom_weights = atom_masses / atom_masses.sum()
        return atom_weights

    def _calculate_com(self, atom_positions: np.ndarray, atom_masses: np.ndarray):
        """
        Calculate the centre of the N atoms to which the iGUESSMD force has been applied, given a set of
        weights defining the weightings of each atom's position.

        :param atom_positions: NumPy array of atom positions with dimensions (N, 3)
        :param atom_masses: NumPy array of atomic masses (AMU) with dimension (N)
        :return: NumPy array containing the position of the centre of mass of the atoms with dimension (3)
        """
        assert self.mobile_atom_indices is not None
        assert np.all(
            atom_positions.shape == np.array((self.mobile_atom_indices.size, 3))
        )
        return calculate_com(atom_positions, atom_masses)

    def _calculate_x_rel_trajectory(self):
        """
        Calculate the trajectory that the COM follows during the iGUESSMD simulation.
        """
        assert (
            self.mobile_atom_indices is not None
            and self.mobile_atom_weights is not None
            and self.reaction_path is not None
        )
        assert not np.array_equal(
            self.iguessmd_sim_mobile_atom_positions,
            np.zeros((self.reaction_path.shape[0], self.mobile_atom_indices.size, 3)),
        )
        self.iguessmd_sim_mobile_centre_positions = np.array(
            [
                calculate_weighted_centre(
                    self.iguessmd_sim_mobile_atom_positions[i],
                    self.mobile_atom_weights,
                )
                for i in range(self.reaction_path.shape[0])
            ]
        )

        if self._coord_type == "abs":
            self.iguessmd_sim_x_rel_positions = (
                self.iguessmd_sim_mobile_centre_positions
            )
            return

        if (
            self._coord_type == "rel"
            and self.ref_atom_indices is not None
            and self.ref_atom_weights is not None
            and self.iguessmd_sim_ref_atom_positions is not None
        ):
            self.iguessmd_sim_ref_centre_positions = np.array(
                [
                    calculate_weighted_centre(
                        self.iguessmd_sim_ref_atom_positions[i],
                        self.ref_atom_weights,
                    )
                    for i in range(self.reaction_path.shape[0])
                ]
            )
            self.iguessmd_sim_x_rel_positions = (
                self.iguessmd_sim_mobile_centre_positions
                - self.iguessmd_sim_ref_centre_positions
            )

    def calculate_cumulative_work_done(self):
        """
        Calculate the cumulative work done by the iGUESSMD force on the
        COM of the atoms with which it interacts over the iGUESSMD simulation.
        """
        self._calculate_x_rel_trajectory()
        self._calculate_iguessmd_forces(self.iguessmd_sim_x_rel_positions)
        self._calculate_work_done()


def iguessmd_com_force(
    parallel_force_constant: float | np.ndarray,
    perpendicular_force_constant: float | np.ndarray,
    uses_pbcs: bool,
    pbc_box_lengths: np.ndarray | None = None,
    coord_type: str = "abs",
):
    """
    Defines a harmonic restraint force for the COM of a group of atoms for performing iGUESSMD.

    :param parallel_force_constant: Force constant of the harmonic restraint to be applied to the
      COM of the group of atoms in the direction tangential to the RC, in units kJ mol-1 nm-2
    :param perpendicular_force_constant: Force constant of the harmonic restraint to be applied to the
      COM of the group of atoms in the direction perpendicular to the RC, in units kJ mol-1 nm-2
    :param uses_pbcs: Bool specifying whether to use periodic boundary conditions for the
      harmonic restraint
    :param pbc_box_lengths: Lengths of the simulation box defining the periodic boundary of the
      simulation
    :return: CustomCentroidBondForce defining the harmonic iGUESSMD force that interacts with
      the COM of the specified atoms.
    """
    assert coord_type in ("abs", "rel")
    n_centres = 1 if coord_type == "abs" else 2

    if uses_pbcs:
        iguessmd_force = CustomCentroidBondForce(
            n_centres, iGUESSMD_FORCE_EXPRESSION_COM_PERIODIC
        )
    else:
        iguessmd_force = CustomCentroidBondForce(
            n_centres, iGUESSMD_FORCE_EXPRESSION_COM_NONPERIODIC
        )

    iguessmd_force.addGlobalParameter(
        f"{iGUESSMD_FORCE_CONSTANT_PARALLEL_PARAMETER_NAME}",
        float(parallel_force_constant),
    )
    iguessmd_force.addGlobalParameter(
        f"{iGUESSMD_FORCE_CONSTANT_PERPENDICULAR_PARAMETER_NAME}",
        float(perpendicular_force_constant),
    )
    if uses_pbcs and pbc_box_lengths is not None:
        iguessmd_force.addGlobalParameter("Lx", pbc_box_lengths[0])
        iguessmd_force.addGlobalParameter("Ly", pbc_box_lengths[1])
        iguessmd_force.addGlobalParameter("Lz", pbc_box_lengths[2])
    if coord_type == "abs":
        iguessmd_force.addGlobalParameter("x2", 0.0)
        iguessmd_force.addGlobalParameter("y2", 0.0)
        iguessmd_force.addGlobalParameter("z2", 0.0)
    iguessmd_force.addPerBondParameter("x0")
    iguessmd_force.addPerBondParameter("y0")
    iguessmd_force.addPerBondParameter("z0")
    iguessmd_force.addPerBondParameter("tx")
    iguessmd_force.addPerBondParameter("ty")
    iguessmd_force.addPerBondParameter("tz")
    iguessmd_force.setUsesPeriodicBoundaryConditions(uses_pbcs)
    iguessmd_force.setForceGroup(31)

    return iguessmd_force


def iguessmd_single_atom_force(
    parallel_force_constant: float | np.ndarray,
    perpendicular_force_constant: float | np.ndarray,
    uses_pbcs: bool,
    pbc_box_lengths: np.ndarray | None = None,
):
    """
    Defines a harmonic restraint force for a single atom for performing iGUESSMD.

    :param parallel_force_constant: Force constant of the harmonic restraint to be applied to the
      specified atom in the direction tangential to the RC, in units kJ mol-1 nm-2
    :param perpendicular_force_constant: Force constant of the harmonic restraint to be applied to the
      specified atom in the direction perpendicular to the RC, in units kJ mol-1 nm-2
    :param uses_pbcs: Bool specifying whether to use periodic boundary conditions for the
      harmonic restraint
    :param pbc_box_lengths: Lengths of the simulation box defining the periodic boundary of the
      simulation
    :return: CustomExternalForce defining the harmonic iGUESSMD force that interacts with the
      specified atom.
    """
    if uses_pbcs:
        iguessmd_force = CustomExternalForce(iGUESSMD_FORCE_EXPRESSION_ATOM_PERIODIC)
    else:
        iguessmd_force = CustomExternalForce(iGUESSMD_FORCE_EXPRESSION_ATOM_NONPERIODIC)
    iguessmd_force.addGlobalParameter(
        f"{iGUESSMD_FORCE_CONSTANT_PARALLEL_PARAMETER_NAME}",
        float(parallel_force_constant),
    )
    iguessmd_force.addGlobalParameter(
        f"{iGUESSMD_FORCE_CONSTANT_PERPENDICULAR_PARAMETER_NAME}",
        float(perpendicular_force_constant),
    )
    if uses_pbcs and pbc_box_lengths is not None:
        iguessmd_force.addGlobalParameter("Lx", pbc_box_lengths[0])
        iguessmd_force.addGlobalParameter("Ly", pbc_box_lengths[1])
        iguessmd_force.addGlobalParameter("Lz", pbc_box_lengths[2])
    iguessmd_force.addPerParticleParameter("x0")
    iguessmd_force.addPerParticleParameter("y0")
    iguessmd_force.addPerParticleParameter("z0")
    iguessmd_force.addPerParticleParameter("tx")
    iguessmd_force.addPerParticleParameter("ty")
    iguessmd_force.addPerParticleParameter("tz")
    iguessmd_force.setForceGroup(31)

    return iguessmd_force
