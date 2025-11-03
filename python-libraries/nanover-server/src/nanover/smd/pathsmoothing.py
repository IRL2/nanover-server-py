import MDAnalysis.units
import numpy as np
import matplotlib.pyplot as plt
import MDAnalysis as mda
import ipywidgets as widgets

from typing import Union, Tuple
from os import PathLike

from scipy.interpolate import splprep, splev

from nanover.mdanalysis import NanoverParser, NanoverReader

try:
    from ipywidgets import interact, interactive, fixed, interact_manual, VBox, HBox
    from IPython import get_ipython
    from IPython.display import display

    _in_notebook = get_ipython().__class__.__name__ in [
        "ZMQInteractiveShell",
        "TerminalInteractiveShell",
    ]
except:
    _in_notebook = False


class PathSmoother:
    """
    A class for smoothing trajectories of atoms or groups of atoms defined
    by user interactions during an iMD simulation.
    """

    @classmethod
    def from_recording(cls, filename: Union[str, PathLike]):
        """
        Create a PathSmoother object from the data stored in a NanoVer recording.

        :param filename: Path to the NanoVer recording.
        """
        ps = cls()
        ps.filename = filename
        ps.retrieve_traj_data()
        (
            ps.atom_positions,
            ps.com_positions,
            ps.atom_indices,
            ps.n_interaction_frames,
        ) = retrieve_uf_com_path(
            ps.all_atom_positions, ps.all_atom_masses, ps.all_user_forces
        )
        return ps

    @classmethod
    def from_path_coordinates(cls, coords: np.ndarray, atom_masses: np.ndarray):
        """
        Create a PathSmoother object from a set of coordinates, atom indices and atom masses that
        directly define the path to be smoothed.

        :param coords: NumPy array with dimensions (m x N x 3), defining the cartesian coordinates of
          the N atoms defining the m frame path to be smoothed, with length units of nanometers
        :param atom_masses: NumPy array with dimensions (N) defining the masses of the atoms for
          which the coordinates have been passed
        """
        ps = cls()
        # Reshape coordinates array if 1-D coordinates given (i.e. if path of single atom given)
        if coords.ndim == 2 and coords.shape[1] == 3:
            ps.atom_positions = np.array([[coords[i]] for i in range(coords.shape[0])])
        elif coords.ndim == 3:
            ps.atom_positions = coords
        ps.atom_masses = atom_masses
        ps.com_positions = np.array(
            [
                calculate_com(ps.atom_positions[i], ps.atom_masses)
                for i in range(ps.atom_positions.shape[0])
            ]
        )
        ps.n_interaction_frames = ps.atom_positions.shape[0]

        return ps

    @staticmethod
    def _make_plots_interactive():
        """
        Function that can be used to make plots interactive in Jupyter or IPython interfaces.
        """
        if _in_notebook:
            get_ipython().run_line_magic("matplotlib", "widget")
        else:
            raise Exception(
                "WARNING: interactive plots are only supported for Jupyter or IPython "
                "interfaces."
            )

    def __init__(self):
        self.filename: Union[str, PathLike] | None = None
        self.com_positions: np.ndarray | None = None
        self.atom_positions: np.ndarray | None = None
        self.atom_indices: np.ndarray | None = None
        self.n_interaction_frames: int | None = None

        self.all_atom_positions: np.ndarray | None = None
        self.all_user_forces: np.ndarray | None = None
        self.all_atom_masses: np.ndarray | None = None

        self.universe: mda.Universe | None = None

        # Plot smoothing properties
        self.fig = None
        self.ax = None
        self.scatter_curve = None
        self.scatter_points = None
        self.smoothing_plot: plt.Figure | None = None
        self.smoothing_parameter: float | None = None
        self.n_points: int | None = None
        self.smoothing_start_index: int | None = None
        self.smoothing_end_index: int | None = None

        self.smoothed_com_trajectory: np.ndarray | None = None
        self.constant_speed_com_trajectory: np.ndarray | None = None
        self.constant_speed_nm_ps: float | None = None
        self.timestep_ps: float | None = None

    def close_interactive_plots(self):
        """
        Function that closes any interactive plots that are currently open within a Jupyter
        or IPython session.
        """
        if _in_notebook:
            plt.close("all")
            print("Existing interactive plots closed.")
            self.fig = None
            self.ax = None
            self.scatter_curve = None
            self.scatter_points = None

    def create_mda_universe(self):
        """
        Create an MDAnalysis universe from the given trajectory file, retaining NanoVer
        units for all quantities.
        """
        assert self.filename is not None
        self.universe = mda.Universe(
            self.filename,
            format=NanoverReader,
            topology_format=NanoverParser,
            convert_units=False,
        )

    def read_mda_universe_data(self):
        """
        Retrieve the relevant data for the trajectory from the MDAnalysis universe.
        """
        assert self.universe is not None

        # Retrieve atom elements and guess masses
        atom_elements = np.array(self.universe.atoms.types)
        # Old code: self.all_atom_masses = mda.topology.guessers.guess_masses(atom_elements)
        self.all_atom_masses = np.array(
            [mda.topology.guessers.get_atom_mass(atom_t) for atom_t in atom_elements]
        )

        # Create empty arrays for atom positions and user forces
        n_atoms = atom_elements.size
        self.all_atom_positions = np.zeros(
            (self.universe.trajectory.n_frames - 1, n_atoms, 3)
        )
        self.all_user_forces = np.zeros(
            (self.universe.trajectory.n_frames - 1, n_atoms, 3)
        )

        # Read trajectory and retrieve atom positions and user forces
        for timestep in self.universe.trajectory:
            # Exclude topology frame where no user forces are applied
            if int(timestep.frame) != 0:
                index = int(timestep.frame) - 1
                self.all_atom_positions[index] = self.universe.atoms.positions
                self.all_user_forces[index] = timestep.data["user_forces"]

    def retrieve_traj_data(self):
        """
        Retrieve data about the trajectory from the given trajectory file.
        """
        self.create_mda_universe()
        self.read_mda_universe_data()

    def plot_com_trajectory(
        self, equal_aspect_ratio: bool = False, cmap: str = "viridis"
    ):
        """
        Plot the trajectory of the COM of the atoms defining the path.

        :param equal_aspect_ratio: A bool defining whether the axes should have equal aspect ratio
        :param cmap: A string defining the Matplotlib colour map to use to plot the trajectory
        """
        assert self.com_positions is not None and self.n_interaction_frames is not None
        self._make_plots_interactive()
        plot_com_trajectory(
            self.com_positions, self.n_interaction_frames, equal_aspect_ratio, cmap
        )

    def plot_atoms_trajectories(
        self, equal_aspect_ratio: bool = False, cmap: str = "viridis"
    ):
        """
        Plot the trajectories of the individual atoms defining the path.

        :param equal_aspect_ratio: A bool defining whether the axes should have equal aspect ratio
        :param cmap: A string defining the Matplotlib colour map to use to plot the trajectory
        """
        assert self.atom_positions is not None and self.n_interaction_frames is not None
        self._make_plots_interactive()
        plot_atom_trajectories(
            self.atom_positions, self.n_interaction_frames, equal_aspect_ratio, cmap
        )

    def create_interactive_smoothing_plot(
        self, equal_aspect_ratio: bool = False, cmap: str = "viridis"
    ):
        """
        Create an interactive plot to smooth the trajectory of the centre of mass of the atoms defining
        the path.

        :param equal_aspect_ratio: A bool defining whether the axes should have equal aspect ratio
        :param cmap: A string defining the Matplotlib colour map to use to plot the trajectory
        """
        self._make_plots_interactive()

        def interactive_smoothing_plot(
            x_pos, y_pos, z_pos, smoothing_value, n_points, start_point, end_point
        ):

            pos_array_size = x_pos.size
            original_u_values = np.linspace(0, 1, pos_array_size)

            interpolated_path, interpolated_u_values = interpolate_path(
                x_pos, y_pos, z_pos, smoothing_value, n_points, start_point, end_point
            )

            # Save current view if plot already exists
            if self.fig is None or self.ax is None:

                # Close the old figure to avoid duplicates
                plt.close(self.fig)

                self.fig = plt.figure(figsize=(8, 8))
                self.ax = plt.axes(111, projection="3d")
                self.scatter_curve = self.ax.scatter3D(
                    interpolated_path[:, 0],
                    interpolated_path[:, 1],
                    interpolated_path[:, 2],
                    c=interpolated_u_values,
                    cmap=cmap,
                    s=1.0,
                )
                self.scatter_points = self.ax.scatter3D(
                    x_pos,
                    y_pos,
                    z_pos,
                    c=original_u_values,
                    cmap=cmap,
                    s=50.0,
                    alpha=0.02,
                )
                self.ax.set_xlabel(r"$x$ / nm")
                self.ax.set_ylabel(r"$y$ / nm")
                self.ax.set_zlabel(r"$z$ / nm")
                if not equal_aspect_ratio:
                    xlim = self.ax.get_xlim3d()
                    ylim = self.ax.get_ylim3d()
                    zlim = self.ax.get_zlim3d()
                    self.ax.set_box_aspect(
                        (xlim[1] - xlim[0], ylim[1] - ylim[0], zlim[1] - zlim[0])
                    )

                cbar = self.fig.colorbar(
                    self.scatter_curve, ax=self.ax, pad=0.1, shrink=0.8
                )
                cbar.set_label("COM trajectory")
                cbar.set_ticks([0, 1])
                cbar.set_ticklabels(
                    [r"$\mathbf{r}_{\rm{start}}$", r"$\mathbf{r}_{\rm{end}}$"]
                )
                plt.show(block=False)

            else:
                self.scatter_curve.remove()
                self.scatter_points.remove()
                self.scatter_curve = self.ax.scatter3D(
                    interpolated_path[:, 0],
                    interpolated_path[:, 1],
                    interpolated_path[:, 2],
                    c=interpolated_u_values,
                    cmap=cmap,
                    s=1.0,
                )
                self.scatter_points = self.ax.scatter3D(
                    x_pos,
                    y_pos,
                    z_pos,
                    c=original_u_values,
                    cmap=cmap,
                    s=50.0,
                    alpha=0.02,
                )

                plt.draw()

            return smoothing_value, n_points, start_point, end_point

        self.smoothing_plot = interactive(
            interactive_smoothing_plot,
            x_pos=fixed(self.com_positions[:, 0]),
            y_pos=fixed(self.com_positions[:, 1]),
            z_pos=fixed(self.com_positions[:, 2]),
            smoothing_value=widgets.FloatSlider(
                min=0.0, max=10.0, step=0.0001, value=0.0, readout_format='.4f'
            ),
            n_points=widgets.IntSlider(min=1000, max=10000, step=100, value=1000),
            start_point=widgets.IntSlider(
                min=0, max=int((self.com_positions[:, 0].size / 2) - 2), step=1, value=0
            ),
            end_point=widgets.IntSlider(
                min=0, max=int((self.com_positions[:, 0].size / 2) - 2), step=1, value=0
            ),
        )

        return self.smoothing_plot

    def save_smoothing_plot_result(self):
        assert self.smoothing_plot is not None
        (
            self.smoothing_parameter,
            self.n_points,
            self.smoothing_start_index,
            self.smoothing_end_index,
        ) = self.smoothing_plot.result
        self.smoothed_com_trajectory, _ = interpolate_path(
            self.com_positions[:, 0],
            self.com_positions[:, 1],
            self.com_positions[:, 2],
            self.smoothing_parameter,
            self.n_points,
            self.smoothing_start_index,
            self.smoothing_end_index,
        )

    def plot_smoothed_com_trajectory(
        self, equal_aspect_ratio: bool = False, cmap: str = "viridis"
    ):
        """
        Plot the smoothed trajectory of the COM of the atoms defining the path.

        :param equal_aspect_ratio: A bool defining whether the axes should have equal aspect ratio
        :param cmap: A string defining the Matplotlib colour map to use to plot the trajectory
        """
        assert self.smoothed_com_trajectory is not None and self.n_points is not None
        self._make_plots_interactive()
        plot_com_trajectory(
            self.smoothed_com_trajectory, self.n_points, equal_aspect_ratio, cmap
        )

    def plot_constant_speed_trajectory(
        self, equal_aspect_ratio: bool = False, cmap: str = "viridis"
    ):
        """
        Plot the constant speed trajectory calculated using :func:`calculate_constant_speed_trajectory`.

        :param equal_aspect_ratio: A bool defining whether the axes should have equal aspect ratio
        :param cmap: A string defining the Matplotlib colour map to use to plot the trajectory
        """
        assert self.constant_speed_com_trajectory is not None
        self._make_plots_interactive()
        n_points = self.constant_speed_com_trajectory.shape[0]
        plot_com_trajectory(
            self.constant_speed_com_trajectory, n_points, equal_aspect_ratio, cmap
        )

    def calculate_constant_speed_trajectory(
        self,
        desired_speed_nm_ps: float,
        timestep_ps: float,
        rel_tolerance: float = 1e-5,
        max_iterations: int = 10,
    ):
        """
        Calculate the coordinates that define a constant speed trajectory along a smoothed
        COM trajectory, given the desired speed in nm ps-1 and timestep of the simulation in ps.

        :param desired_speed_nm_ps: The desired speed in nm ps-1
        :param timestep_ps: The timestep in ps
        :param rel_tolerance: The relative tolerance in the errors of the speed and distance between points
        :param max_iterations: The maximum number of iterations for the iterative refinement procedure
        """

        assert self.smoothed_com_trajectory is not None
        path_length = calculate_trajectory_length(self.smoothed_com_trajectory)
        n_points_required = int(
            np.ceil(path_length / (desired_speed_nm_ps * timestep_ps))
        )  # + 1
        refined_path, _ = interpolate_path(
            self.smoothed_com_trajectory[:, 0],
            self.smoothed_com_trajectory[:, 1],
            self.smoothed_com_trajectory[:, 2],
            0.0,
            n_points=n_points_required,
        )

        print(
            f"Starting iterative path refinement for trajectory of length {path_length} nm with desired speed {desired_speed_nm_ps} nm...\n"
        )

        # Iterate until convergence criteria achieved
        for i in range(max_iterations):
            print(f"\n--------------\nIteration {i+1}\n--------------\n")
            n_points_required = int(
                np.ceil(path_length / (desired_speed_nm_ps * timestep_ps))
            )  # + 1
            print(
                f"Number of points required to achieve desired speed along path: {n_points_required}\n"
            )
            refined_path, _ = interpolate_path(
                refined_path[:, 0],
                refined_path[:, 1],
                refined_path[:, 2],
                0.0,
                n_points_required,
            )
            print("Path interpolation completed\n")
            path_length = calculate_trajectory_length(refined_path)
            speed_along_path_nm_ps = path_length / (n_points_required * timestep_ps)
            print(f"Total path length: {path_length} nm")
            print(f"Speed along path: {speed_along_path_nm_ps} nm ps-1")
            relative_speed_error = (
                abs(speed_along_path_nm_ps - desired_speed_nm_ps)
            ) / desired_speed_nm_ps
            distances = np.linalg.norm(np.diff(refined_path, axis=0), axis=1)
            maximum_distance_error = (
                np.max(distances) - np.min(distances)
            ) / np.average(distances)
            if (
                maximum_distance_error < rel_tolerance
                and relative_speed_error < rel_tolerance
            ):
                print(
                    "***\nRelative tolerance achieved, exiting iterative refinement...\n***\n\n"
                )
                print(f"Final path length: {path_length} nm")
                print(f"Final speed: {speed_along_path_nm_ps} nm ps-1")
                self.constant_speed_com_trajectory = refined_path
                self.constant_speed_nm_ps = speed_along_path_nm_ps
                self.timestep_ps = timestep_ps
                break
            elif i == max_iterations - 1:
                print(
                    "***\nFailed to achieve relative tolerance, aborting iterative refinement...\n***\n\n"
                )
                break

        return None

    def save_constant_speed_trajectory_data(self, output_filepath: PathLike[str]):
        """
        Save the constant speed trajectory data as numpy arrays to an output file
        whose path is defined by the user.
        """
        assert output_filepath is not None
        assert (
            self.constant_speed_com_trajectory is not None
            and self.constant_speed_nm_ps is not None
            and self.timestep_ps is not None
        )

        print(f"Saving constant speed trajectory data to {output_filepath}")

        with open(output_filepath, "wb") as outfile:
            np.save(outfile, self.constant_speed_com_trajectory)
            np.save(outfile, self.constant_speed_nm_ps)
            np.save(outfile, self.timestep_ps)
            if self.atom_indices is not None:
                np.save(outfile, self.atom_indices)
            else:
                print(
                    "Warning: atom indices not provided and hence not saved to output file."
                )


def get_uf_atoms_and_frames(user_forces: np.ndarray) -> np.ndarray:
    """
    Takes a numpy array of user forces and returns a numpy array containing the indices of
    the atoms to which the user forces were applied and the indices of the frames during which
    the user forces were applied.

    :param user_forces: A NumPy array of user forces with dimensions (n_frames, n_particles, 3).
    :return indices_and_frames: An array of vectors containing the atom index (first element) and
     frame index (second element) encoding to which atom and when the user forces were applied.
    """

    return np.unique(np.transpose(np.nonzero(user_forces))[:, 0:2], axis=0)


def get_uf_first_indices_atom_frame(
    indices_and_frames: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Takes a numpy array of atom indices to which user forces are applied and the indices of the
    frames during which they are applied, and returns the index of the first frame in which the
    first user force was applied and the indices of the atoms to which it was applied.

    :param indices_and_frames: The NumPy array output by :func:`get_uf_atoms_and_frames`.
    :return Tuple[first_frame_index, atom_indices]: A tuple of the frame index in which the first
     user force was applied, and an array of indices of the atoms to which it was applied.
    """

    first_user_force_frame_index = indices_and_frames[0, 0]
    user_force_atom_indices = indices_and_frames[
        np.where(indices_and_frames[:, 0] == first_user_force_frame_index)[0]
    ][:, 1]

    return first_user_force_frame_index, user_force_atom_indices


def separate_atom_and_frame_indices(
    indices_and_frames: np.ndarray,
) -> Tuple[np.ndarray, list]:
    """
    Takes the array output by :func:`get_uf_atoms_and_frames` and returns separate arrays of atom
    indices and frame indices defining the frames in which those forces were applied.

    :param indices_and_frames: Array output by :func:`get_uf_atoms_and_frames`
    :return Tuple[atom_indices, frame_indices]: Tuple containing an array of atom indices and an
     array of frame indices.
    """

    indices = np.unique(indices_and_frames[:, 1])
    frames_for_indices = []
    for index in indices:
        frames_for_indices.append(
            indices_and_frames[np.where(indices_and_frames[:, 1] == index)][:, 0]
        )
    return indices, frames_for_indices


def check_continuous_user_force(
    atom_indices: np.ndarray,
    frames_for_indices: list,
    user_force_atom_indices: np.ndarray,
) -> None:
    """
    A function that checks that the user force for which path smoothing will occur was continuously
    applied to the system.

    :param atom_indices: A NumPy array of the indices of atoms to which user forces were applied.
    :param frames_for_indices: A list of NumPy arrays defining the frames in which the user forces were applied to
     the atoms in the corresponding atom_indices array.
    :param user_force_atom_indices: The indices of the atoms to which the first user force was applied.
    """

    for user_force_atom_index in user_force_atom_indices:
        assert np.all(
            np.diff(
                frames_for_indices[
                    np.where(atom_indices == user_force_atom_index)[0][0]
                ]
            )[:]
            == 1
        )


def get_uf_final_index_frame(
    atom_indices: np.ndarray,
    frames_for_indices: list,
    user_force_atom_indices: np.ndarray,
) -> np.int64:
    """
    Takes the array of indices to which user forces were applied in the simulation, the corresponding list of arrays
    of frames to which those user forces were applied, and the indices of the atoms to which the first user force was
    first applied, and returns the index of the final frame in which the first user force was applied.

    :param atom_indices: A NumPy array of the indices of atoms to which user forces were applied.
    :param frames_for_indices: A list of NumPy arrays defining the frames in which the user forces were applied to
     the atoms in the corresponding atom_indices array.
    :param user_force_atom_indices: The indices of the atoms to which the first user force was applied.
    :return final_user_force_frame_index: The index of the final frame in which the first user force was applied.
    """

    final_indices = np.zeros(user_force_atom_indices.size)
    for i in range(user_force_atom_indices.size):
        final_indices[i] = frames_for_indices[
            np.where(atom_indices == user_force_atom_indices[i])[0][0]
        ][-1]

    # Check that all final indices agree
    (final_indices == final_indices[0]).all()

    return np.int64(final_indices[0])


def get_indices_defining_path(
    user_forces: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray, np.int64]:
    """
    Takes a NumPy array of user forces and returns the indices of the atoms to which the first user force was applied,
    and the indices of the first and final frames in which that user force was applied.

    :param user_forces: A NumPy array of user forces with dimensions (n_frames, n_particles, 3).
    :return Tuple[atom_indices, first_frame_index, last_frame_index]: A tuple containing the indices of the atoms to
     which the first user force was applied, and the indices of the first and final frames in which that user force
     was applied.
    """

    indices_and_frames = get_uf_atoms_and_frames(user_forces)
    first_frame_index, atom_indices = get_uf_first_indices_atom_frame(
        indices_and_frames
    )
    all_atom_indices, frames_for_indices = separate_atom_and_frame_indices(
        indices_and_frames
    )
    check_continuous_user_force(all_atom_indices, frames_for_indices, atom_indices)
    final_frame_index = get_uf_final_index_frame(
        all_atom_indices, frames_for_indices, atom_indices
    )

    return atom_indices, first_frame_index, final_frame_index


def retrieve_uf_com_path(
    all_atom_traj_positions: np.ndarray,
    all_atom_masses: np.ndarray,
    all_user_forces: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.int64]:
    """
    Takes a numpy array of the positions of all atoms across the full trajectory, the masses of the atoms and the full
    array of user forces applied during the trajectory, and returns the path of the COM of the atoms to which the first user force
    was applied across the frames in which that force was applied.

    :param all_atom_traj_positions: A NumPy array defining the positions of all of the atoms across the trajectory.
    :param all_atom_masses: A NumPy array defining the masses of the atoms.
    :param all_user_forces: A NumPy array defining all user forces applied to the atoms throughout the trajectory.
    :return Tuple[atom_path, atom_index, n_frames]: A tuple containing a NumPy array of the positions of the COM of
     the atoms to which the user first applied force that defines the path the COM of the atoms took during the
     application of the force, the indices of the atoms and the number of frames in which that force was applied.
    """

    atom_indices, first_frame_index, final_frame_index = get_indices_defining_path(
        all_user_forces
    )

    atom_masses = all_atom_masses[atom_indices]

    atom_positions = all_atom_traj_positions[
        first_frame_index:final_frame_index, atom_indices, :
    ]
    n_frames = atom_positions.shape[0]

    com_positions = np.array(
        [calculate_com(atom_positions[i], atom_masses) for i in range(n_frames)]
    )

    return atom_positions, com_positions, atom_indices, n_frames


def retrieve_uf_atoms_path(
    all_atom_traj_positions: np.ndarray, all_user_forces: np.ndarray
) -> Tuple[np.ndarray, np.ndarray, np.int64]:
    """
    Takes a numpy array of the positions of all atoms across the full trajectory and the full array of user
    forces applied during the trajectory, and returns the paths of the atoms to which the first user force
    was applied across the frames in which that force was applied.

    :param all_atom_traj_positions: A NumPy array defining the positions of all of the atoms across the trajectory.
    :param all_user_forces: A NumPy array defining all user forces applied to the atoms throughout the trajectory.
    :return Tuple[atoms_path, atom_index, n_frames]: A tuple containing a NumPy array of the positions of the atoms
     to which the user first applied force that define the paths the atoms took during the application of the force,
     the indices of the atoms and the number of frames in which that force was applied.
    """

    atoms_indices, first_frame_index, final_frame_index = get_indices_defining_path(
        all_user_forces
    )

    atoms_positions = all_atom_traj_positions[
        first_frame_index:final_frame_index, atoms_indices, :
    ]
    n_frames = atoms_positions.shape[0]

    return atoms_positions, atoms_indices, n_frames


def calculate_com(atom_positions: np.ndarray, atom_masses: np.ndarray) -> np.ndarray:
    """
    Calculate the COM of a group of atoms given the positions of the atoms and their masses

    :param atom_positions: A n_atoms x 3 NumPy array of atom positions for a given frame
    :param atom_masses: A n_atoms NumPy array of atom masses
    :return np.ndarray: A NumPy array defining the position of the COM of the atoms
    """

    return np.sum(
        np.multiply(np.transpose(atom_positions), atom_masses), axis=1
    ) / np.sum(atom_masses)


def plot_com_trajectory(
    atom_positions: np.ndarray,
    n_frames: int,
    equal_aspect_ratio: bool = False,
    cmap: str = "viridis",
) -> None:
    """
    Function that takes the trajectory of an atom as a NumPy array and the number of frames of the trajectory,
    and create a 3-D plot of the trajectory, coloured using the Viridis colour scheme.

    :param atom_positions: A NumPy array of atom positions defining the trajectory of the atom
    :param n_frames: An integer defining the number of frames of the trajectory
    :param equal_aspect_ratio: A bool defining whether the axes should have equal aspect ratio
    :param cmap: A string defining the Matplotlib colour map to use to plot the trajectory
    """
    assert atom_positions.shape[0] == n_frames
    ax = plt.axes(projection="3d")
    ax.scatter3D(
        atom_positions[:, 0],
        atom_positions[:, 1],
        atom_positions[:, 2],
        c=np.linspace(0, 1, n_frames),
        cmap=cmap,
        s=1.0,
    )
    ax.set_xlabel(r"$x$ / nm")
    ax.set_ylabel(r"$y$ / nm")
    ax.set_zlabel(r"$z$ / nm")
    if not equal_aspect_ratio:
        xlim = ax.get_xlim3d()
        ylim = ax.get_ylim3d()
        zlim = ax.get_zlim3d()
        ax.set_box_aspect((xlim[1] - xlim[0], ylim[1] - ylim[0], zlim[1] - zlim[0]))
    plt.show()


def plot_atom_trajectories(
    atoms_positions: np.ndarray,
    n_frames: int,
    equal_aspect_ratio: bool = False,
    cmap: str = "viridis",
) -> None:
    """
    Function that takes the trajectory of an atom as a NumPy array and the number of frames of the trajectory,
    and create a 3-D plot of the trajectory, coloured using the Viridis colour scheme.

    :param atoms_positions: A NumPy array of atom positions defining the trajectory of the atom
    :param n_frames: An integer defining the number of frames of the trajectory
    :param equal_aspect_ratio: A bool defining whether the axes should have equal aspect ratio
    :param cmap: A string defining the Matplotlib colour map to use to plot the trajectory
    """
    assert atoms_positions.shape[0] == n_frames
    ax = plt.axes(projection="3d")
    for i in range(atoms_positions.shape[1]):
        ax.scatter3D(
            atoms_positions[:, i, 0],
            atoms_positions[:, i, 1],
            atoms_positions[:, i, 2],
            c=np.linspace(0, 1, n_frames),
            cmap=cmap,
            s=1.0,
        )
    ax.set_xlabel(r"$x$ / nm")
    ax.set_ylabel(r"$y$ / nm")
    ax.set_zlabel(r"$z$ / nm")
    if not equal_aspect_ratio:
        xlim = ax.get_xlim3d()
        ylim = ax.get_ylim3d()
        zlim = ax.get_zlim3d()
        ax.set_box_aspect((xlim[1] - xlim[0], ylim[1] - ylim[0], zlim[1] - zlim[0]))
    plt.show()


def interpolate_path(
    x_pos: np.ndarray,
    y_pos: np.ndarray,
    z_pos: np.ndarray,
    smoothing_parameter: float,
    n_points: int,
    start_index: int = 0,
    end_index: int = 0,
):
    """
    Function that interpolates the path defined by the arrays of x, y and z coordinates passed to
    it and returns an n_point trajectory of coordinates that are approximately equally spaced along
    the interpolated path.

    :param x_pos: A 1-D NumPy array of x coordinates
    :param y_pos: A 1-D NumPy array of y coordinates
    :param z_pos: A 1-D NumPy array of z coordinates
    :param smoothing_parameter: A float defining the smoothing parameter to be used in the interpolation
    :param n_points: An integer defining the number of interpolation points
    :param start_index: An integer defining the index of the coordinates at which the interpolation
      should begin.
    :param end_index: An integer defining the index of the coordinates at which the interpolation
      should end.
    """
    pos_array_size = x_pos.size
    tck, u = splprep(
        [
            x_pos[start_index : pos_array_size - end_index],
            y_pos[start_index : pos_array_size - end_index],
            z_pos[start_index : pos_array_size - end_index],
        ],
        s=smoothing_parameter,
    )
    u_fine = np.linspace(0, 1, n_points)
    x_fine, y_fine, z_fine = splev(u_fine, tck)
    return np.transpose(np.squeeze([x_fine, y_fine, z_fine])), u_fine


def calculate_trajectory_length(path: np.ndarray):
    """
    Calculate the total length of an m frame trajectory defined by an mx3 array of
    coordinates in nanometers.

    :param path: A NumPy array of coordinates with dimensions m x 3
    """

    return np.sum(np.linalg.norm(np.diff(path, axis=0), axis=1))
