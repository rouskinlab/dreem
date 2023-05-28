"""
DREEM RNAstructure module.

Wrapper around RNAstructure from the Mathews Lab at U of Rochester:
https://rna.urmc.rochester.edu/RNAstructure.html
"""

from logging import getLogger
from pathlib import Path

from ..core.rna import RnaProfile
from ..core.shell import run_cmd

logger = getLogger(__name__)


def fold(rna: RnaProfile, *,
         out_dir: Path, temp_dir: Path, save_temp: bool,
         dms_quantile: float,
         rerun: bool):
    """ Run the 'Fold' program of RNAstructure. """
    ct_file = rna.get_ct(out_dir)
    if rerun or not ct_file.is_file():
        cmd = ["Fold"]
        if dms_quantile > 0.:
            # Write the DMS reactivities file for the RNA.
            cmd.extend(["--DMS", rna.to_dms(temp_dir, dms_quantile)])
        # Write a temporary FASTA file for the RNA.
        cmd.append(fasta := rna.to_fasta(temp_dir))
        try:
            # Get the path of the output CT file.
            cmd.append(ct_file)
            # Run the command.
            run_cmd(cmd, check_is_after=[ct_file])
        finally:
            if not save_temp:
                # Delete the temporary files.
                fasta.unlink(missing_ok=True)
    else:
        logger.warning(f"File already exists: {ct_file}")
    return ct_file
