from os import PathLike
from pathlib import Path
from typing import Any, Optional, Sequence, Tuple

import numpy as np

from openmm.app import Simulation, StateDataReporter
from openmm.unit import nanometer

from nanover.core import AppServer

from .converter import openmm_to_frame_data
from . import serializer
from .imd import (
    create_imd_force,
    add_imd_force_to_system,
    ImdForceManager,
    NON_IMD_FORCES_GROUP_MASK,
)
from .thermo import compute_instantaneous_temperature, compute_dof
from nanover.imd.imd_force import calculate_contribution_to_work

from nanover.trajectory.frame_wrapper import FrameData

class OpenMMSimulation:
    """
    A wrapper for OpenMM simulations to run inside the OmniRunner.

    The following attributes can be configured after construction:
    - :attr:`frame_interval`: Number of simulation steps to advance between frames.
    - :attr:`include_velocities`: Include particle velocities in frames.
    - :attr:`include_forces`: Include particle forces in frames.
    - :attr:`platform_name`: Name of OpenMM platform to use when loading the system from XML.
    """

    @classmethod
    def from_simulation(cls, simulation: Simulation, *, name: str | None = None):
        """
        Construct this from an existing OpenMM simulation.

        :param simulation: An existing OpenMM Simulation
        :param name: An optional name for the simulation instead of default
        """
        sim = cls(name)
        sim.simulation = simulation
        sim.imd_force = add_imd_force_to_system(simulation.system)
        sim.simulation.context.reinitialize(preserveState=True)
        sim.determine_pbcs()

        sim.checkpoint = sim.simulation.context.createCheckpoint()

        return sim

    @classmethod
    def from_xml_path(cls, path: PathLike[str], *, name: str | None = None):
        """
        Construct this from an existing NanoVer OpenMM XML file at a given path.

        :param path: Path of the NanoVer OpenMM XML file
        :param name: An optional name for the simulation instead of filename
        """
        if name is None:
            name = Path(path).stem
        sim = cls(name)
        sim.xml_path = path
        return sim

    def __init__(self, name: str | None = None):
        self.name = name or "Unnamed OpenMM Simulation"

        self.xml_path: PathLike[str] | None = None
        self.app_server: AppServer | None = None

        self.frame_interval = 5
        """Number of simulation steps to advance between frames."""
        self.include_velocities = False
        """Include particle velocities in frames."""
        self.include_forces = False
        """Include particle forces in frames."""
        self.platform_name: str | None = None
        """Name of OpenMM platform to use at the time the system is loaded from XML."""
        self.use_pbc_wrapping: bool | None = None
        """Provide atom positions wrapped according to PBC such that each molecule has a center of mass within the
        primary periodic box."""
        self.pbc_vectors: np.ndarray | None = None
        """Array of vectors defining the periodic box used by the simulation (if PBCs are employed)."""

        self.imd_force = create_imd_force()
        self.simulation: Simulation | None = None
        self.checkpoint: Any | None = None
        self.verbose_reporter: StateDataReporter | None = None

        self.imd_force_manager: ImdForceManager | None = None

        self.work_done: float = 0.0
        self._work_done_intermediate: float = 0.0
        self._prev_imd_forces: np.ndarray | None = None
        self._prev_imd_indices: np.ndarray | None = None

        self._dof: int | None = None

    def load(self):
        """
        Load and set up the simulation if it isn't done already.
        """
        if self.xml_path is None or self.simulation is not None:
            return

        with open(self.xml_path) as infile:
            self.imd_force = create_imd_force()
            self.simulation = serializer.deserialize_simulation(
                infile, imd_force=self.imd_force, platform_name=self.platform_name
            )

        self.determine_pbcs()
        self.checkpoint = self.simulation.context.createCheckpoint()

    def reset(self, app_server: AppServer):
        """
        Reset the simulation to its initial conditions, reset IMD interactions, and reset frame stream to begin with
        topology and continue.

        :param app_server: The app server hosting the frame publisher and imd state
        """
        assert self.simulation is not None and self.checkpoint is not None

        self.app_server = app_server
        self.imd_force_manager = ImdForceManager(
            self.app_server.imd,
            self.imd_force,
            self.pbc_vectors if self.use_pbc_wrapping else None,
        )

        self._dof = compute_dof(self.simulation.system)

        # reset imd and work
        self.work_done = 0.0
        self._work_done_intermediate = 0.0
        self._prev_imd_forces = None
        self._prev_imd_indices = None

        # reload initial state and cleanup forces
        self.simulation.context.reinitialize()
        self.simulation.context.loadCheckpoint(self.checkpoint)

        # send the initial topology frame
        frame_data = self.make_topology_frame()
        self.app_server.frame_publisher.send_clear()
        self.app_server.frame_publisher.send_frame(frame_data)

        # verbose reporter
        if (
            self.verbose_reporter is not None
            and self.verbose_reporter not in self.simulation.reporters
        ):
            self.simulation.reporters.append(self.verbose_reporter)

    def determine_pbcs(self):
        """
        Determine whether the simulation uses periodic boundary conditions and if it does,
        retrieve the periodic box vectors in nanometers.
        """
        assert self.simulation is not None

        if self.use_pbc_wrapping is False:
            return

        self.use_pbc_wrapping = self.simulation.system.usesPeriodicBoundaryConditions()
        if self.use_pbc_wrapping:
            self.pbc_vectors = np.array(
                [
                    vector.value_in_unit(nanometer)
                    for vector in self.simulation.system.getDefaultPeriodicBoxVectors()
                ]
            )

    def advance_by_seconds(self, dt: float):
        """
        Advance playback time by some seconds, and advance the simulation to the next frame output.

        :param dt: Time to advance playback by in seconds (ignored)
        """
        self.advance_to_next_report()

    def advance_by_one_step(self):
        """
        Advance the simulation to the next point a frame should be reported, and send that frame.
        """
        self.advance_to_next_report()

    def advance_to_next_report(self):
        """
        Step the simulation to the next point a frame should be reported, and send that frame.
        """
        assert (
            self.simulation is not None
            and self.imd_force_manager is not None
            and self.app_server is not None
        )

        # determine step count for next frame
        steps_to_next_frame = (
            self.frame_interval - self.simulation.currentStep % self.frame_interval
        )

        # advance the simulation
        self.simulation.step(steps_to_next_frame)

        # fetch positions early, for updating imd
        state = self.simulation.context.getState(
            getPositions=True,
            enforcePeriodicBox=self.use_pbc_wrapping or False,
        )
        positions = state.getPositions(asNumpy=True)

        # Calculate on-step contribution to work
        if self._prev_imd_forces is not None:
            affected_atom_positions = positions[self._prev_imd_indices]
            self._work_done_intermediate += calculate_contribution_to_work(
                self._prev_imd_forces, affected_atom_positions
            )

        # update imd forces and energies
        self.imd_force_manager.update_interactions(self.simulation, positions)

        # generate the next frame with the existing (still valid) positions
        frame_data = self.make_regular_frame(positions)

        # Update work done in frame data
        self.work_done = self._work_done_intermediate
        frame_data.user_work_done = self.work_done

        # Calculate previous-step contribution to work for the next time step
        # (negative contribution, so subtract from the total work done)
        if frame_data.user_forces_sparse is not None:
            affected_atom_positions = positions[frame_data.user_forces_index]
            self._work_done_intermediate -= calculate_contribution_to_work(
                frame_data.user_forces_sparse, affected_atom_positions
            )

        # send the next frame
        self.app_server.frame_publisher.send_frame(frame_data)

        # Update previous step forces (saving them in their sparse form)
        self._prev_imd_forces = frame_data.user_forces_sparse
        self._prev_imd_indices = frame_data.user_forces_index

    def make_topology_frame(self):
        """
        Make a NanoVer FrameData corresponding to the current particle positions and topology of the simulation.
        """
        assert self.simulation is not None

        state = self.simulation.context.getState(
            getPositions=True,
            getEnergy=True,
            enforcePeriodicBox=self.use_pbc_wrapping or False,
        )
        topology = self.simulation.topology
        frame_data = openmm_to_frame_data(state=state, topology=topology)
        return frame_data

    def make_regular_frame(self, positions: np.ndarray | None = None):
        """
        Make a NanoVer FrameData corresponding to the current state of the simulation.

        :param positions: Optionally provided particle positions to save fetching them again.
        """
        assert (
            self.simulation is not None
            and self.imd_force_manager is not None
            and self._dof is not None
        )

        # fetch omm state
        state = self.simulation.context.getState(
            getPositions=positions is None,
            getForces=self.include_forces,
            getVelocities=self.include_velocities,
            getEnergy=True,
            enforcePeriodicBox=self.use_pbc_wrapping or False,
            groups=NON_IMD_FORCES_GROUP_MASK,
        )

        # generate frame based on basic omm state
        frame_data = openmm_to_frame_data(
            state=state,
            topology=None,
            include_positions=positions is None,
            include_velocities=self.include_velocities,
            include_forces=self.include_forces,
            state_excludes_imd=True,
        )

        # Assume that the KE is always available, which is true for this case
        frame_data.system_temperature = compute_instantaneous_temperature(
            self.simulation, frame_data.kinetic_energy, self._dof
        )

        # add any provided positions
        if positions is not None:
            frame_data.particle_positions = positions.astype(np.float32)

        # add imd force information
        self.imd_force_manager.add_to_frame_data(frame_data)

        return frame_data

class HBondOpenMMSimulation(OpenMMSimulation):
    def __init__(self,*args , **kwargs):
        super().__init__(*args, **kwargs)
        self._cfg_h_range: Optional[Tuple[int, int]] = None
        self._cfg_a_range: Optional[Tuple[int, int]] = None
        self._cfg_h_indices: Optional[np.ndarray] = None  
        self._cfg_a_indices: Optional[np.ndarray] = None
        self._particle_elements: Optional[np.ndarray] = None 
        self._hydrogen_indices: Optional[np.ndarray] = None  
        self._nof_indices: Optional[np.ndarray] = None  
        self._last_hb_pairs: set[tuple[int, int]] = set()

    def make_regular_frame(self, positions: np.ndarray | None = None):
        """
        Make a NanoVer FrameData corresponding to the current state of the simulation.

        :param positions: Optionally provided particle positions to save fetching them again.
        """
        assert (
            self.simulation is not None
            and self.imd_force_manager is not None
            and self._dof is not None
        )

        # fetch omm state
        state = self.simulation.context.getState(
            getPositions=positions is None,
            getForces=self.include_forces,
            getVelocities=self.include_velocities,
            getEnergy=True,
            enforcePeriodicBox=self.use_pbc_wrapping or False,
            groups=NON_IMD_FORCES_GROUP_MASK,
        )

        # generate frame based on basic omm state
        frame_data = openmm_to_frame_data(
            state=state,
            topology=self.simulation.topology,
            include_positions=positions is None,
            include_velocities=self.include_velocities,
            include_forces=self.include_forces,
            state_excludes_imd=True,
        )

        # Assume that the KE is always available, which is true for this case
        frame_data.system_temperature = compute_instantaneous_temperature(
            self.simulation, frame_data.kinetic_energy, self._dof
        )

        # add any provided positions
        if positions is not None:
            frame_data.particle_positions = positions.astype(np.float32)

        # Add hydrogen bond calculation here
        self.updateHydrogenbond(frame_data, cutoff=0.3)

        # add imd force information
        self.imd_force_manager.add_to_frame_data(frame_data)

        return frame_data
    
    def set_hbond_selection(self,
        *,
        h_indices: Optional[Sequence[int]] = None,
        acceptor_indices: Optional[Sequence[int]] = None,
        h_index_range: Optional[Tuple[int, int]] = None,
        a_index_range: Optional[Tuple[int, int]] = None,
        ) -> None:
        def _as_intp(arr: Optional[Sequence[int]]) -> Optional[np.ndarray]:
            if arr is None:
                return None
            return np.asarray(arr, dtype=np.intp)
        self._cfg_h_indices = _as_intp(h_indices)
        self._cfg_a_indices = _as_intp(acceptor_indices)
        self._cfg_h_range   = h_index_range
        self._cfg_a_range   = a_index_range
        self._hydrogen_indices = None
        self._nof_indices = None

    def _ensure_element_and_indices(self, N_hint: Optional[int] = None) -> bool:
        if self._hydrogen_indices is not None and self._nof_indices is not None:
            return (self._hydrogen_indices.size > 0 and self._nof_indices.size > 0)
        #check if we have the hydrogen and acceptor indices
        has_h = (self._cfg_h_indices is not None) or (self._cfg_h_range is not None)
        has_a = (self._cfg_a_indices is not None) or (self._cfg_a_range is not None)
        if not (has_h and has_a):
            return False  
        
        #if we have ranges, we need to get the element types
        need_elements = (self._cfg_h_range is not None) or (self._cfg_a_range is not None)
        if need_elements and self._particle_elements is None:
            if self.simulation is None or self.simulation.topology is None:
                return False
            topo = self.simulation.topology
            elems = [a.element.atomic_number if a.element is not None else 0 for a in topo.atoms()]
            self._particle_elements = np.asarray(elems, dtype=np.uint8)

        N = N_hint if N_hint is not None else (
            self._particle_elements.size if self._particle_elements is not None else None
        )

        #hydrogen 
        if self._cfg_h_indices is not None:
            H = self._cfg_h_indices
            if N is not None:
                H = H[(H >= 0) & (H < N)]
            self._hydrogen_indices = H
        elif self._cfg_h_range is not None:
            if self._particle_elements is None:
                return False
            lo, hi = self._cfg_h_range
            lo = max(0, min(lo, self._particle_elements.size))
            hi = max(lo, min(hi, self._particle_elements.size))
            H_local = np.where(self._particle_elements[lo:hi] == 1)[0]  # H=1
            self._hydrogen_indices = (H_local + lo).astype(np.intp, copy=False)
        else:
            return False

        # acceptor (N, O, F)
        if self._cfg_a_indices is not None:
            A = self._cfg_a_indices
            if N is not None:
                A = A[(A >= 0) & (A < N)]
            self._nof_indices = A
        elif self._cfg_a_range is not None:
            if self._particle_elements is None:
                return False
            lo, hi = self._cfg_a_range
            lo = max(0, min(lo, self._particle_elements.size))
            hi = max(lo, min(hi, self._particle_elements.size))
            A_local = np.where(np.isin(self._particle_elements[lo:hi], [7, 8, 9]))[0]  # N/O/F
            self._nof_indices = (A_local + lo).astype(np.intp, copy=False)
        else:
            return False

        return (self._hydrogen_indices.size > 0 and self._nof_indices.size > 0)


    def updateHydrogenbond(self, data: FrameData, cutoff: float):
        particle_positions = getattr(data, "particle_positions", None)
        if particle_positions is None:
            return
        positions = np.asarray(particle_positions, dtype=np.float32)

        if not self._ensure_element_and_indices(N_hint=positions.shape[0]):
            return

        H = self._hydrogen_indices
        A = self._nof_indices

        h_pos = positions[H] #shape (nH, 3)  
        a_pos = positions[A] #shape (nA, 3)
        d = np.linalg.norm(h_pos[:, None, :] - a_pos[None, :, :], axis=2)  
        jmin = d.argmin(axis=1)
        dmin = d[np.arange(H.size), jmin]
        ok = dmin <= cutoff

        hb_pairs = set()
        for idx in np.where(ok)[0]:
            hi = int(H[idx]); aj = int(A[jmin[idx]])
            pair = (hi, aj) if hi < aj else (aj, hi)
            hb_pairs.add(pair)

        existing_pairs = np.asarray(data.bond_pairs, dtype=np.uint32).reshape(-1, 2)
        existing_set = {(int(i), int(j)) for i, j in existing_pairs}
        if not hasattr(self, "_last_hb_pairs"):
            self._last_hb_pairs = set()
        keep_set = existing_set - self._last_hb_pairs
        keep_mask  = [(int(i), int(j)) in keep_set for i, j in existing_pairs]
        kept_pairs = existing_pairs[keep_mask]
        new_pairs = sorted(hb_pairs - keep_set)
        if new_pairs:
            new_pairs = np.asarray(new_pairs, dtype=np.uint32).reshape(-1, 2) 
        else:
            new_pairs = np.empty((0, 2), dtype=np.uint32)  
        data.bond_pairs = np.vstack([kept_pairs, new_pairs]).astype(np.uint32)
        self._last_hb_pairs = set(map(tuple, new_pairs.tolist()))
    