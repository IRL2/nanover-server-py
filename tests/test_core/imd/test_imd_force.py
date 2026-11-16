import hypothesis
import numpy as np
import pytest
from hypothesis import strategies, given
from math import exp
from nanover.imd.imd_force import (
    get_center_of_mass_subset,
    calculate_spring_force,
    calculate_gaussian_force,
    apply_single_interaction_force,
    calculate_imd_force,
    calculate_constant_force,
    InvalidInteractionError, INTERACTION_METHOD_MAP, ForceCalculator,
)
from nanover.imd.particle_interaction import ParticleInteraction
from nanover.testing.strategies import vec3s

# precomputed results of gaussian force.
EXP_1 = exp(-1 / 2)
EXP_3 = exp(-3 / 2)
UNIT = np.array([1, 1, 1]) / np.linalg.norm([1, 1, 1])


@pytest.fixture
def particles():
    num_particles = 50
    positions = np.array([[i, i, i] for i in range(num_particles)])
    masses = np.array([i + 1 for i in range(num_particles)])
    return positions, masses


@pytest.fixture
def single_interaction():
    position = (0, 0, 0)
    index = 1
    return ParticleInteraction(
        position=position,
        particles=[index],
    )


@pytest.fixture
def single_interaction_multiple_atoms():
    position = (0, 0, 0)
    return ParticleInteraction(
        position=position,
        particles=[1, 2, 3],
    )


def test_multiple_interactions(particles):
    """
    Tests multiple concurrent interactions.

    Ensures that equidistant interactions on particles [0,1] and particles [1,2] results in zero force on particle 1,
    and the same (but opposite) forces on atoms 0 and 2.
    """

    positions, masses = particles
    interaction = ParticleInteraction(
        position=[0.5, 0.5, 0.5],
        particles=[0, 1],
    )
    interaction_2 = ParticleInteraction(
        position=[1.5, 1.5, 1.5],
        particles=[1, 2],
    )
    # set masses of atoms 0 and 2 to be the same, so things cancel out nicely.
    masses[2] = masses[0]
    single_forces = np.zeros((len(positions), 3))

    single_energy = apply_single_interaction_force(
        positions, masses, interaction, single_forces
    )

    energy, forces = calculate_imd_force(
        positions, masses, [interaction, interaction_2]
    )
    expected_energy = 2 * single_energy
    expected_forces = np.zeros((len(positions), 3))
    expected_forces[0, :] = single_forces[0, :]
    expected_forces[1, :] = 0
    # the uneven masses in this calculation mean both atoms 0 and 2 are pulled towards atom 1.
    expected_forces[2, :] = -single_forces[0, :]

    assert np.allclose(energy, expected_energy)
    assert np.allclose(forces, expected_forces)


@pytest.mark.parametrize("particle_count", (1, 2))
def test_interaction_invalid_particle_index(particles, particle_count):
    """
    Test that attempting to calculate iMD forces for interactions with out of bounds particles raises InvalidInteractionError.
    """
    positions, masses = particles
    single_forces = np.zeros((len(positions), 3))
    indexes = [len(positions) + i for i in range(particle_count)]

    interaction = ParticleInteraction(particles=indexes)

    with pytest.raises(InvalidInteractionError):
        apply_single_interaction_force(positions, masses, interaction, single_forces)

    with pytest.raises(InvalidInteractionError):
        calculate_imd_force(positions, masses, [interaction])


@pytest.mark.parametrize("scale", [np.nan, np.inf, -np.inf])
def test_interaction_force_invalid_scale(particles, single_interaction, scale):
    with pytest.raises(ValueError):
        single_interaction.scale = scale


@pytest.mark.parametrize("scale", [-1.0, 0, 100])
def test_interaction_force_single(particles, single_interaction, scale):
    """
    Tests that the interaction force calculation gives the expected result on a single atom, at a particular position,
    with varying scale.
    """
    positions, masses = particles
    forces = np.zeros((len(positions), 3))
    expected_forces = np.zeros((len(positions), 3))
    single_interaction.scale = scale
    energy = apply_single_interaction_force(
        positions, masses, single_interaction, forces
    )

    diff = positions[1, :] - single_interaction.position
    expected_energy = (1 - EXP_3) * scale
    expected_forces[1, :] = np.array(
        [
            - diff
            * EXP_3
            * scale
            * (
                masses[single_interaction.particles[0]]
                / np.sum(masses[single_interaction.particles[0]])
            )
        ]
    )

    assert np.allclose(energy, expected_energy, equal_nan=True)
    assert np.allclose(forces, expected_forces, equal_nan=True)


@pytest.mark.parametrize("max_force", [np.nan, -np.inf, -1])
def test_invalid_max_force(single_interaction, max_force):
    with pytest.raises(ValueError):
        single_interaction.max_force = max_force


# TODO: does it make any sense to test NaN, infinite, and negative masses?
@pytest.mark.parametrize("mass", [-1.0, 100, np.nan, np.inf, -np.inf])
def test_interaction_force_mass(
    particles,
    single_interaction,
    mass,
):
    """
    tests that the interaction force calculation gives the expected result on a single atom, at a particular position,
    with varying mass.
    """
    positions, masses = particles
    forces = np.zeros((len(positions), 3))
    expected_forces = np.zeros((len(positions), 3))
    masses = np.array([mass] * len(masses))
    energy = apply_single_interaction_force(
        positions, masses, single_interaction, forces
    )

    expected_energy = 1 - EXP_3
    diff = positions[1, :] - single_interaction.position
    expected_forces[1, :] = np.array([- diff * EXP_3 * (mass / mass)])

    assert np.allclose(energy, expected_energy, equal_nan=True)
    assert np.allclose(forces, expected_forces, equal_nan=True)


def test_interaction_force_zero_mass_singleatom(particles, single_interaction):
    positions, masses = particles
    forces = np.zeros((len(positions), 3))
    masses = np.array([0.0] * len(masses))

    energy = apply_single_interaction_force(
        positions, masses, single_interaction, forces
    )
    assert energy == pytest.approx(0)


def test_interaction_ignores_massless(
    particles, single_interaction_multiple_atoms
):
    """
    Tests that inclusion or exclusion of massless particles in the interaction does not change the resulting energy and
    forces.
    """
    positions, masses = particles

    # make first particle massless
    masses[0] = 0

    # without first particle
    single_interaction_multiple_atoms.particles = [1, 2, 3]
    forces_A = np.zeros((len(positions), 3))
    energy_A = apply_single_interaction_force(
        positions, masses, single_interaction_multiple_atoms, forces_A
    )

    # with first particle
    single_interaction_multiple_atoms.particles = [0, 1, 2, 3]
    forces_B = np.zeros((len(positions), 3))
    energy_B = apply_single_interaction_force(
        positions, masses, single_interaction_multiple_atoms, forces_B
    )

    assert energy_A == energy_B
    assert np.allclose(forces_A, forces_B)


def test_interaction_force_zero_mass_multiatom(
    particles, single_interaction_multiple_atoms
):
    positions, masses = particles
    forces = np.zeros((len(positions), 3))
    masses = np.array([0.0] * len(masses))

    with pytest.raises(ZeroDivisionError):
        apply_single_interaction_force(
            positions, masses, single_interaction_multiple_atoms, forces
        )


@pytest.mark.parametrize(
    "position, selection, selection_masses",
    [
        ([1, 1, 1], [0, 1], [1, 2]),
        ([2, 2, 2], [0, 1], [1, 2]),
        ([0, 0, 0], [0, 1], [1, 2]),
        ([0, 0, 0], [0, 1, 49], [1, 2, 0]),
        ([0, 0, 0], [0, 1, 49], [1, 2, 10]),
        ([-5, -5, -5], [0, 1, 49], [1, 2, 10]),
        ([np.nan, np.nan, np.nan], [0, 1], [1, 2]),
    ],
)
def test_interaction_force_com(particles, position, selection, selection_masses):
    """
    tests that the interaction force gives the correct result when acting on a group of atoms.
    """
    position = np.array(position)
    selection = np.array(selection)
    interaction = ParticleInteraction(
        position=position,
        particles=selection,
    )
    positions, masses = particles
    # set non uniform masses based on parameterisation
    for index, mass in zip(selection, selection_masses):
        masses[index] = mass
    forces = np.zeros((len(positions), 3))

    # perform the full calculation to generate expected result.
    com = get_center_of_mass_subset(positions, masses, selection)
    diff = com - interaction.position
    dist_sqr = np.dot(diff, diff)
    exponential = exp(-dist_sqr / 2)
    expected_energy = 1 - exponential
    expected_forces = np.zeros((len(positions), 3))
    selection_mass = np.sum(masses[selection])
    for index in selection:
        expected_forces[index, :] = (
            -diff * masses[index] / selection_mass * exponential
        )

    energy = apply_single_interaction_force(positions, masses, interaction, forces)
    assert np.allclose(energy, expected_energy, equal_nan=True)
    assert np.allclose(forces, expected_forces, equal_nan=True)


@pytest.mark.parametrize(
    "position, selection, selection_masses",
    [
        ([1, 1, 1], [0, 1], [1, 2]),
        ([2, 2, 2], [0, 1], [1, 2]),
        ([0, 0, 0], [0, 1], [1, 2]),
        ([0, 0, 0], [0, 1, 49], [1, 2, 0]),
        ([0, 0, 0], [0, 1, 49], [1, 2, 10]),
        ([-5, -5, -5], [0, 1, 49], [1, 2, 10]),
        ([np.nan, np.nan, np.nan], [0, 1], [1, 2]),
    ],
)
def test_interaction_force_no_mass_weighting(
    particles, position, selection, selection_masses
):
    """
    tests that the interaction force gives the correct result when acting on a group of atoms.
    """
    position = np.array(position)
    selection = np.array(selection)
    interaction = ParticleInteraction(
        position=position,
        particles=selection,
        mass_weighted=False,
    )
    positions, masses = particles
    # Set non uniform masses based on parameterisation
    for index, mass in zip(selection, selection_masses):
        masses[index] = mass
    forces = np.zeros((len(positions), 3))

    # Perform explicit calculation to find expected energy
    com = get_center_of_mass_subset(positions, masses, selection)
    diff = com - interaction.position
    dist_sqr = np.dot(diff, diff)
    exponential = np.exp(-dist_sqr / 2)
    expected_energy = 1 - exponential

    # Calculate normalised weights for user forces (and energies)
    weights = (masses[selection] != 0).astype(int)
    weights = weights / np.sum(weights)

    # Calculate expected forces
    expected_forces = np.zeros((len(positions), 3))
    expected_forces[selection, :] = - exponential * np.outer(weights, diff)

    # Retrieve and check energy and forces
    energy = apply_single_interaction_force(positions, masses, interaction, forces)
    assert np.allclose(energy, expected_energy, equal_nan=True)
    assert np.allclose(forces, expected_forces, equal_nan=True)


def test_interaction_force_unknown_type(particles, single_interaction):
    single_interaction.interaction_type = "unknown_type"
    positions, masses = particles
    forces = np.zeros((len(positions), 3))

    with pytest.raises(KeyError):
        apply_single_interaction_force(positions, masses, single_interaction, forces)


def test_get_com_all(particles):
    positions, masses = particles
    subset = [i for i in range(len(positions))]

    com = get_center_of_mass_subset(positions, masses, subset)

    expected_com = np.sum(positions * masses[:, None], axis=0) / masses.sum()

    assert np.allclose(com, expected_com)


def test_get_com_subset(particles):
    positions, masses = particles
    subset = [i for i in range(0, len(positions), 2)]

    com = get_center_of_mass_subset(positions, masses, subset)

    expected_com = (
        np.sum(positions[subset] * masses[subset, None], axis=0) / masses[subset].sum()
    )

    assert np.allclose(com, expected_com)


@strategies.composite
def random_periodic_box_lengths(draw):
    # Generate random polar coordinates and convert them to euclidean
    # coordinates to get a periodic box.
    # box length has to be nonzero.
    length = strategies.floats(
        min_value=0.01, max_value=100,
        allow_nan=False, allow_infinity=False
    )

    return np.array([draw(length) for _ in range(3)])


@strategies.composite
def random_positions_pbc(draw, particle_count=10):
    periodic_box_lengths = draw(random_periodic_box_lengths())

    # pick two random points in lowest quadrant of the box.
    # TODO positions at or very near zero cause problems, as the wrap can flip between 0 and box length.
    lengths = np.array(
        [
            strategies.floats(
                min_value=0.005,
                max_value=box_length * 0.5,
                allow_nan=False,
                allow_infinity=False,
            )
            for box_length in periodic_box_lengths
        ]
    )

    random_masses = strategies.floats(
        min_value=0.01, max_value=100,
        allow_nan=False, allow_infinity=False
    )
    masses = np.array([draw(random_masses) for _ in range(particle_count)])

    positions = np.zeros((particle_count, 3))
    for i in range(particle_count):
        positions[i] = np.array([draw(coord) for coord in lengths])

    # generate random integer values to multiply positions by, putting them in different images.
    images = strategies.integers(min_value=-100, max_value=100)
    image_multiples = np.array([draw(images) for _ in range(3 * particle_count)])
    image_multiples.reshape((particle_count, 3))

    # move points to new random positions around the periodic box.
    positions_periodic = np.zeros((particle_count, 3))
    for i in range(particle_count):
        positions_periodic[i] = positions[i] + image_multiples[i] * periodic_box_lengths

    return positions, masses, positions_periodic, periodic_box_lengths


@strategies.composite
def random_positions_pbc_subset(draw, particle_count=10):
    positions, masses, positions_periodic, periodic_box_lengths = draw(random_positions_pbc(particle_count))
    particles = draw(strategies.lists(strategies.integers(min_value=0, max_value=particle_count-1), unique=True))
    return particles, positions, masses, positions_periodic, periodic_box_lengths


@given(random_positions_pbc())
def test_get_com_subset_pbc(positions_pbc):
    """
    Tests that the center of mass calculation works when using periodic boundary conditions.
    """
    positions, masses, positions_periodic, periodic_box_lengths = positions_pbc
    subset = [i for i in range(0, len(positions), 2)]

    expected_com = (
        np.sum(positions[subset] * masses[subset, None], axis=0) / masses[subset].sum()
    )

    com = get_center_of_mass_subset(
        positions_periodic,
        masses,
        subset,
        periodic_box_lengths=periodic_box_lengths,
    )

    assert np.allclose(com, expected_com)


def test_get_com_single():
    position = np.array([[1, 0, 0]])
    mass = np.array([20])
    com = get_center_of_mass_subset(position, mass, subset=[0])
    assert np.allclose(com, position)


def test_get_com_different_array_lengths(particles):
    positions, mass = particles
    # make masses array not match positions in length.
    mass = np.array([1])
    subset = range(positions.shape[0])
    with pytest.raises(IndexError):
        get_center_of_mass_subset(positions, mass, subset)


@pytest.mark.parametrize(
    "position, interaction_position, expected_energy, expected_force",
    [
        ([1, 0, 0], [0, 0, 0], 1 - EXP_1, [-EXP_1, 0, 0]),
        ([0, 0, 0], [1, 0, 0], 1 - EXP_1, [EXP_1, 0, 0]),
        ([1, 3, 0], [1, 2, 0], 1 - EXP_1, [0, -EXP_1, 0]),
        ([1, 3, 3], [1, 3, 2], 1 - EXP_1, [0, 0, -EXP_1]),
        (UNIT, [0, 0, 0], 1 - EXP_1, np.multiply(UNIT, [-EXP_1, -EXP_1, -EXP_1])),
        ([1, 2, 3], [1, 2, 3], 1 - 1, [0, 0, 0]),
        ([1, 1, 1], [0, 0, 0], 1 - EXP_3, [-EXP_3] * 3),
        ([1, 0, 0], [1, 0, 0], 1 - 1, [0, 0, 0]),
        ([-1, -1, -1], [0, 0, 0], 1 - EXP_3, [EXP_3] * 3),
    ],
)
def test_gaussian_force(
    position, interaction_position, expected_energy, expected_force
):
    # tests gaussian force for various hand evaluated values.
    energy, force = calculate_gaussian_force(
        np.array(position), np.array(interaction_position)
    )
    assert np.allclose(energy, expected_energy, equal_nan=True)
    assert np.allclose(force, expected_force, equal_nan=True)


@pytest.mark.parametrize(
    "position, interaction, expected_energy, expected_force",
    [
        ([1, 0, 0], [0, 0, 0], 1, [-2, 0, 0]),
        ([0, 0, 0], [1, 0, 0], 1, [2, 0, 0]),
        ([1, 3, 0], [1, 2, 0], 1, [0, -2, 0]),
        ([1, 3, 3], [1, 3, 2], 1, [0, 0, -2]),
        (UNIT, [0, 0, 0], 1, np.multiply(UNIT, [-2, -2, -2])),
        ([1, 1, 1], [0, 0, 0], 3, [-2, -2, -2]),
        ([1, 2, 3], [1, 2, 3], 0, [0, 0, 0]),
        ([-1, -1, -1], [0, 0, 0], 3, [2, 2, 2]),
    ],
)
def test_spring_force(position, interaction, expected_energy, expected_force):
    energy, force = calculate_spring_force(np.array(position), np.array(interaction))
    assert np.allclose(energy, expected_energy, equal_nan=True)
    assert np.allclose(force, expected_force, equal_nan=True)


CONSTANT_TESTS = [
    (
        position,
        interaction,
        np.linalg.norm(np.subtract(interaction, position)),
        np.subtract(interaction, position)
        / np.linalg.norm(np.subtract(interaction, position)),
    )
    for (position, interaction) in [
        ([1, 0, 0], [0, 0, 0]),
        ([0, 0, 0], [1, 0, 0]),
        ([1, 3, 0], [1, 2, 0]),
        ([1, 3, 3], [1, 3, 2]),
        (UNIT, [0, 0, 0]),
        ([1, 1, 1], [0, 0, 0]),
        ([-1, -1, -1], [0, 0, 0]),
    ]
]


@pytest.mark.parametrize(
    "position, interaction, expected_energy, expected_force",
    CONSTANT_TESTS,
)
def test_constant_force(position, interaction, expected_energy, expected_force):
    energy, force = calculate_constant_force(np.array(position), np.array(interaction))
    assert np.allclose(energy, expected_energy, equal_nan=True)
    assert np.allclose(force, expected_force, equal_nan=True)

INTERACTION_SCALES = strategies.floats(width=32, min_value=0, allow_infinity=False, allow_nan=False)
INTERACTION_MAX_FORCES = strategies.floats(width=32, min_value=0, allow_nan=False)
INTERACTION_TYPES = strategies.sampled_from(list(INTERACTION_METHOD_MAP.keys()))

FORCE_CALCULATORS = list(INTERACTION_METHOD_MAP.values())

@hypothesis.given(
    position=vec3s(),
    calculator=strategies.sampled_from(FORCE_CALCULATORS),
    periodic_box_lengths=random_periodic_box_lengths(),
    force_magnitude_limit=INTERACTION_MAX_FORCES,
)
def test_overlap_no_force(
    position,
    calculator: ForceCalculator,
    periodic_box_lengths,
    force_magnitude_limit: float,
):
    """Test that overlapping interaction and position results in no force for all interaction types."""
    energy, force = calculator(
        particle_position=np.array(position),
        interaction_position=np.array(position),
        periodic_box_lengths=periodic_box_lengths,
        force_magnitude_limit=force_magnitude_limit,
    )
    assert energy == 0
    assert np.allclose(force, 0)


@hypothesis.given(
    system=random_positions_pbc_subset(),
    interaction_position=vec3s(),
    interaction_type=INTERACTION_TYPES,
    scale=INTERACTION_SCALES,
    max_force=INTERACTION_MAX_FORCES,
)
def test_interaction_force_max_force(
    system,
    interaction_position,
    interaction_type: str,
    scale: float,
    max_force: float,
):
    """Test that total force applied never exceeds max force for all interaction types."""

    particles, positions, masses, positions_periodic, periodic_box_lengths = system

    forces = np.zeros((len(positions), 3))

    interaction = ParticleInteraction(
        position=interaction_position,
        scale=scale,
        max_force=max_force,
        particles=particles,
        interaction_type=interaction_type,
    )

    _ = apply_single_interaction_force(
        positions=positions,
        masses=masses,
        interaction=interaction,
        forces=forces,
        periodic_box_lengths=periodic_box_lengths,
    )

    force_magnitude = np.sum(np.linalg.norm(forces, axis=1))
    assert force_magnitude <= max_force or np.isclose(force_magnitude, max_force)


@hypothesis.given(
    particle_position=vec3s(),
    interaction_position=vec3s(),
    calculator=strategies.sampled_from(FORCE_CALCULATORS),
    periodic_box_lengths=random_periodic_box_lengths(),
    force_magnitude_limit=INTERACTION_MAX_FORCES,
)
def test_force_magnitude_limit(
    particle_position,
    interaction_position,
    calculator: ForceCalculator,
    periodic_box_lengths,
    force_magnitude_limit: float,
):
    """Test that the force magnitude limit is never exceeded for all interaction types."""
    _, force = calculator(
        particle_position=np.array(particle_position),
        interaction_position=np.array(interaction_position),
        periodic_box_lengths=periodic_box_lengths,
        force_magnitude_limit=force_magnitude_limit,
    )
    force_magnitude = np.sum(np.linalg.norm(force))
    assert force_magnitude <= force_magnitude_limit or np.isclose(force_magnitude, force_magnitude_limit)
