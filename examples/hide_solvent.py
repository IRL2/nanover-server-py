from nanover.app import NanoverImdClient

SOLVENT_RESIDUE_NAME = "HOH"

with NanoverImdClient.autoconnect() as client:
    client.subscribe_to_frames()
    client.subscribe_multiplayer()
    first_frame = client.wait_until_first_frame(check_interval=0.5, timeout=10)

    print(f"Attempting to hide residue {SOLVENT_RESIDUE_NAME} (residues in frame: {', '.join(set(first_frame.residue_names))})")

    particles = (
        particle_index 
        for particle_index, residue_index 
        in enumerate(first_frame.particle_residues) 
        if first_frame.residue_names[residue_index] == SOLVENT_RESIDUE_NAME
    )

    solvent = client.create_selection("solvent", particles)
    with solvent.modify():
        solvent.hide = True
        solvent.interaction_method = "none"
