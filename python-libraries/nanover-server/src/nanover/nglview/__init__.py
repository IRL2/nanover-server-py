"""
Module providing an NGLView python client for use in Jupyter notebooks or IPython.
"""

from nanover.jupyter import NGLClient as NGLClient

import warnings

warnings.warn(
    "use `from nanover.jupyter import NGLClient`",
    DeprecationWarning,
    stacklevel=2,
)
