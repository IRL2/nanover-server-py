import numpy as np
import msgpack
from typing import Iterable
from nanover.mdanalysis.converter import ELEMENT_INDEX

def write_webtraj(selection, io):
    universe = selection.universe
    
    topology = {
        "elements": pack_array("B", (ELEMENT_INDEX[e] for e in selection.elements)),
        "bonds": pack_array("I", selection.bonds.to_indices().flat),
    }
    
    positions = []
    for t in universe.trajectory:
        positions.append(pack_array("f", (c * .1 for c in selection.positions.flat)))
    
    data = {
        "topology": topology,
        "positions": positions,
    }

    io.write(msgpack.packb(data))

def pack_array(dtype: str, values: Iterable):
    return np.fromiter(values, dtype=dtype).tobytes()
