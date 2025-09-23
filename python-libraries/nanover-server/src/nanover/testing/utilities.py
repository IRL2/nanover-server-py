from typing import Any
import numpy as np


def simplify_numpy(dict_frame: dict[str, Any]):
    """
    Convert in-place all numpy arrays into simple lists so pytest can check equality naively, and return the dict for
    convenience.
    :param dict_frame:
    """
    for key, value in dict_frame.items():
        if isinstance(value, np.ndarray):
            dict_frame[key] = list(value)
    return dict_frame
