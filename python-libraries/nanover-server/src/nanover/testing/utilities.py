from typing import Any
import numpy as np


def simplify_numpy(dict_frame: dict[str, Any]):
    """
    Copy dict with numpy arrays converted into simple lists.

    Useful for pytest equality checks and some conversions.
    :param dict_frame:
    """
    simplified = {}
    for key, value in dict_frame.items():
        if isinstance(value, np.ndarray):
            simplified[key] = value.reshape(-1).tolist()
        else:
            simplified[key] = value
    return simplified
