"""
DREEM RNAstructure module.

Functions to run RNAstructure from the Mathews Lab at U of Rochester:
https://rna.urmc.rochester.edu/RNAstructure.html
"""

from pathlib import Path

from ..core import path
from ..core.rna import RnaProfile
from ..core.shell import run_cmd


def fold(out_dir: Path, temp_dir: Path, *,
         rna: RnaProfile,
         dms: bool,
         percentile: float):
    cmd = ["Fold"]
    if dms:
        # Write the DMS reactivities file for the RNA.
        cmd.extend(["--DMS", rna.to_dms(temp_dir, percentile)])
    # Write a temporary FASTA file for the RNA.
    cmd.append(rna.to_fasta(temp_dir))
    # Get the path of the output CT file.
    cmd.append(ct := rna.get_ct(out_dir))
    # Run the command.
    run_cmd(cmd, verify_outputs=[ct])
    return ct
