"""
DREEM RNAstructure module.

Functions to run RNAstructure from the Mathews Lab at U of Rochester:
https://rna.urmc.rochester.edu/RNAstructure.html
"""

from abc import ABC, abstractmethod
from pathlib import Path

from ..core import path
from ..core.rna import RnaProfile
from ..core.shell import run_cmd


class Command(ABC):
    """ Base class for RNAstructure commands. """
    def __init__(self, out_dir: Path, temp_dir: Path):
        self.out_dir = out_dir
        self.temp_dir = temp_dir

    @abstractmethod
    def run(self):
        """ Run the RNAstructure command. """
        return


class Fold(Command):
    def __init__(self, out_dir: Path, temp_dir: Path, *,
                 rna: RnaProfile,
                 max_dms: float = 0.0):
        super().__init__(out_dir, temp_dir)
        self.rna = rna
        self.max_dms = max_dms

    def run(self):
        cmd = ["Fold"]
        if self.max_dms > 0.:
            # Write the DMS reactivities file for the RNA.
            cmd.extend(["--DMS", self.rna.to_dms(self.temp_dir, self.max_dms)])
        # Write a temporary FASTA file for the RNA.
        cmd.append(self.rna.to_fasta(self.temp_dir))
        # Get the path of the output CT file.
        cmd.append(ct := self.rna.get_ct(self.temp_dir))
        # Run the command.
        run_cmd(cmd, verify_outputs=[ct])
        return ct
