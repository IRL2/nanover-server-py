import warnings

from nanover.client import NanoverImdClient as NanoverImdClient

warnings.warn(
    "use 'from nanover.client import NanoverImdClient' instead",
    DeprecationWarning,
    stacklevel=2,
)
