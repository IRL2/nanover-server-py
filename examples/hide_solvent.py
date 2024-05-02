"""
Example of using selections to hide solvent by residue name.
"""

from nanover.app import NanoverImdClient
from nanover.mdanalysis import frame_data_to_mdanalysis
from nanover.trajectory import FrameData

SOLVENT_RESIDUE_NAME = "HOH"


def get_selection_indices(frame: FrameData, query: str):
    universe = frame_data_to_mdanalysis(frame)
    atoms = universe.select_atoms(query)
    indices = map(int, atoms.indices)
    return indices


with NanoverImdClient.autoconnect() as client:
    client.subscribe_to_frames()
    client.subscribe_multiplayer()
    first_frame = client.wait_until_first_frame(check_interval=0.5, timeout=10)

    print(f"Attempting to hide residue {SOLVENT_RESIDUE_NAME} (residues in frame: {', '.join(set(first_frame.residue_names))})")

    solvent_indices = get_selection_indices(client.frame, f"resname {SOLVENT_RESIDUE_NAME}")
    solvent_selection = client.create_selection("solvent", solvent_indices)
    with solvent_selection.modify():
        solvent_selection.hide = True
        solvent_selection.interaction_method = "none"
