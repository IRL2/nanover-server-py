"""
Example of using selections to hide solvent by residue name.
"""

from nanover.app import NanoverImdClient
from nanover.mdanalysis import frame_data_to_mdanalysis
from nanover.trajectory import FrameData

SOLVENT_RESIDUE_NAME = "HOH"


def get_selection_indices(frame: FrameData, query: str):
    # use frame's topology to construct an mdanalysis universe
    universe = frame_data_to_mdanalysis(frame)
    # query mdanalysis universe for desired atoms
    atoms = universe.select_atoms(query)
    # convert to integer atom indices
    indices = map(int, atoms.indices)
    return indices


with NanoverImdClient.autoconnect() as client:
    # we need frames so we can query topology and multiplayer to set shared selections
    client.subscribe_to_frames()
    client.subscribe_multiplayer()

    # wait for an initial frame in which topology will be available
    first_frame = client.wait_until_first_frame(check_interval=0.5, timeout=10)

    print(f"Attempting to hide residue {SOLVENT_RESIDUE_NAME} (residues in frame: {', '.join(set(first_frame.residue_names))})")

    # get atom indicies matching an mdanalysis selection for the particular residue name
    solvent_indices = get_selection_indices(first_frame, f"resname {SOLVENT_RESIDUE_NAME}")

    # create the selection with the desired particles hidden and non-interactable
    with client.create_selection("solvent").modify() as selection:
        selection.set_particles(solvent_indices)
        selection.hide = True
        selection.interaction_method = "none"
