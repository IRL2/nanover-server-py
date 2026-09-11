from io import BytesIO
from zipfile import ZipFile

import pytest
from openmm.unit import nanometer
from simulation_utils import (
    basic_simulation_bundle,
    build_basic_simulation,
    empty_imd_force,
)

from nanover.app import NanoverImdApplication
from nanover.openmm import OpenMMSimulation
from nanover.openmm.bundles import unbundle_openmm_simulation, bundle_openmm_simulation
from nanover.openmm.imd import create_imd_force


def remove_bundle_file(bundle: BytesIO, file_to_remove: str) -> BytesIO:
    next_bundle = BytesIO()

    with (
        ZipFile(bundle, mode="r") as prev,
        ZipFile(next_bundle, mode="w") as next,
    ):
        for file in prev.filelist:
            if file.filename != file_to_remove:
                next.writestr(file.filename, prev.open(file.filename).read())

    return next_bundle


@pytest.mark.parametrize("imd_force", (None, create_imd_force()))
def test_bundle_runs(basic_simulation_bundle, imd_force):
    """Test that a simulation bundle can be unbundled and run."""
    simulation = unbundle_openmm_simulation(
        basic_simulation_bundle, imd_force=imd_force
    )
    assert simulation.currentStep == 0
    simulation.step(5)
    assert simulation.currentStep == 5


@pytest.mark.parametrize(
    "file_missing", ("topology.pdbx", "system.xml", "integrator.xml")
)
def test_missing_file(basic_simulation_bundle, file_missing):
    """Test that the unbundling raises an exception when a required file is missing."""
    broken_bundle = remove_bundle_file(basic_simulation_bundle, file_missing)
    with pytest.raises(ValueError):
        unbundle_openmm_simulation(broken_bundle)


def test_imd_force(basic_simulation_bundle, empty_imd_force):
    """Test that when unbundling a simulation while passing an imd force, the force is added to the system."""
    simulation = unbundle_openmm_simulation(
        basic_simulation_bundle, imd_force=empty_imd_force
    )
    # The provided force is appended as the last step, so will be the last listed force
    force_obtained = simulation.system.getForces()[-1]
    force_added = empty_imd_force
    # The forces are the same if by modifying one we also modify the other.
    force_added.setParticleParameters(0, 0, (1.0, 2.0, 3.0))
    parameters = force_obtained.getParticleParameters(0)
    assert parameters == [0, (1.0, 2.0, 3.0)]


@pytest.mark.parametrize("platform", ("Reference", "CPU"))
def test_platform(basic_simulation_bundle, platform):
    """Test that we can choose the platform when deserialising a simulation."""
    simulation = unbundle_openmm_simulation(
        basic_simulation_bundle, platform_name=platform
    )
    effective_platform_name = simulation.context.getPlatform().getName()
    assert effective_platform_name == platform


def test_serializer_pbc():
    """Test that the periodic boundary conditions are correctly bundled."""
    omm_sim = build_basic_simulation()
    UNIT_SIMULATION_BOX_VECTORS = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    omm_sim.context.setPeriodicBoxVectors(*UNIT_SIMULATION_BOX_VECTORS)

    def out_of_bounds(coord):
        return not 0 <= coord <= 1

    def get_sim_position_coords(sim):
        for position in sim.make_regular_frame().particle_positions:
            yield from position

    with BytesIO() as bundle_pbc:
        bundle_openmm_simulation(omm_sim, outfile=bundle_pbc, pbc_wrapping=True)
        sim_pbc_obj = unbundle_openmm_simulation(bundle_pbc)

    openmm_sim = OpenMMSimulation.from_simulation(sim_pbc_obj)
    openmm_sim.load()
    with NanoverImdApplication.basic_server(port=0) as app_server:
        openmm_sim.reset(app_server)

    assert not any(
        out_of_bounds(coord) for coord in get_sim_position_coords(openmm_sim)
    )

    with BytesIO() as bundle_no_pbc:
        bundle_openmm_simulation(omm_sim, outfile=bundle_no_pbc, pbc_wrapping=False)
        sim_no_pbc_obj = unbundle_openmm_simulation(bundle_no_pbc)

    openmm_sim = OpenMMSimulation.from_simulation(sim_no_pbc_obj)
    openmm_sim.load()
    with NanoverImdApplication.basic_server(port=0) as app_server:
        openmm_sim.reset(app_server)

    assert any(out_of_bounds(coord) for coord in get_sim_position_coords(openmm_sim))
