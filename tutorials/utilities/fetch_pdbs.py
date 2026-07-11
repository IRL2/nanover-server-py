"""
This functionality will exist in the next version of MDAnalysis.
"""

try:
    import pooch
except ImportError:
    print("Requires pooch library (`pip install pooch`)")

import MDAnalysis as mda
from pathlib import Path


def fetch_pdb_universes(*pdb_ids, to_guess=("types", "masses")):
    filenames = {pdb_id: f"{pdb_id}.pdb.gz" for pdb_id in pdb_ids}
    registry_dictionary = {filename: None for filename in filenames.values()}

    downloader = pooch.create(
        path="downloaded_pdbs",
        base_url="https://files.wwpdb.org/download/",
        registry=registry_dictionary,
    )

    return {
        pdb_id: mda.Universe(
            Path(downloader.fetch(fname=file_name)),
            to_guess=to_guess,
        )
        for pdb_id, file_name in filenames.items()
    }
