import warnings

warnings.warn(
    "use `from nanover.app import OmniRunner`",
    DeprecationWarning,
    stacklevel=2,
)

from nanover.app import OmniRunner
