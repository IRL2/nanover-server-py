"""
Keys used in NanoVer simulation frames.
See: https://irl2.github.io/nanover-docs/concepts/applications.html#the-trajectory-application
"""

FRAME_INDEX = "frame.index"

BOX_VECTORS = "system.box.vectors"
SIMULATION_TIME = "system.simulation.time"
SIMULATION_COUNTER = "system.simulation.counter"
SIMULATION_EXCEPTION = "system.simulation.exception"

BOND_PAIRS = "bond.pairs"
BOND_ORDERS = "bond.orders"

PARTICLE_POSITIONS = "particle.positions"
PARTICLE_VELOCITIES = "particle.velocities"
PARTICLE_FORCES = "particle.forces"
PARTICLE_FORCES_SYSTEM = "particle.forces.system"
PARTICLE_ELEMENTS = "particle.elements"
PARTICLE_NAMES = "particle.names"
PARTICLE_RESIDUES = (
    "particle.residues"  # Index of the residue each particle belongs to.
)
PARTICLE_COUNT = "particle.count"

RESIDUE_NAMES = "residue.names"
RESIDUE_IDS = "residue.ids"  # Index of the chain each residue belongs to.
RESIDUE_CHAINS = "residue.chains"
RESIDUE_COUNT = "residue.count"

CHAIN_NAMES = "chain.names"
CHAIN_COUNT = "chain.count"

KINETIC_ENERGY = "energy.kinetic"
POTENTIAL_ENERGY = "energy.potential"
USER_ENERGY = "energy.user.total"

USER_FORCES_SPARSE = "forces.user.sparse"
USER_FORCES_INDEX = "forces.user.index"
USER_WORK_DONE = "forces.user.work_done"

SYSTEM_TEMPERATURE = "system.temperature"

SERVER_TIMESTAMP = "server.timestamp"

PLAY_COMMAND = "playback/play"
RESET_COMMAND = "playback/reset"
STEP_COMMAND = "playback/step"
STEP_BACK_COMMAND = "playback/step_back"
PAUSE_COMMAND = "playback/pause"
LIST_COMMAND = "playback/list"
LOAD_COMMAND = "playback/load"
NEXT_COMMAND = "playback/next"
