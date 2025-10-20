import warnings

warnings.warn(
    "use `from nanover.ase import ASESimulation`",
    DeprecationWarning,
    stacklevel=2,
)

from nanover.ase import ASESimulation
