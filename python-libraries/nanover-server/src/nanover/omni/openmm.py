import warnings

warnings.warn(
    "use `from nanover.openmm import OpenMMSimulation`",
    DeprecationWarning,
    stacklevel=2,
)

from nanover.openmm import OpenMMSimulation
