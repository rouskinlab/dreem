import os
from typing import List, Tuple

import pandas as pd
import subprocess
import warnings

from dreem.util.util import DNA, run_cmd
from dreem.util.cli import OUT_DIR, LIBRARY, COORDS, PRIMERS, FILL, PARALLEL
from dreem.util.path import BAM_EXT
from dreem.vector.mprofile import VectorWriterSpawner
from dreem.util.files_sanity import check_library


def add_coords_from_library(library_path: str,
                            coords: List[Tuple[str, int, int]],
                            fasta: str,
                            out_dir: str):
    library = check_library(pd.read_csv(library_path), fasta, out_dir)
    for row in library.index:
        construct = os.path.splitext(str(library.loc[row, "construct"]))[0]
        first = int(library.loc[row, "section_start"])
        last = int(library.loc[row, "section_end"])
        coord = (construct, first, last)
        if coord not in coords:
            coords.append(coord)


def encode_primers(primers: List[Tuple[str, str, str]]):
    return [(ref, DNA(fwd.encode()), DNA(rev.encode()))
            for ref, fwd, rev in primers]


def run(fasta: str, bam_dirs: List[str], out_dir: str = OUT_DIR,
        library: str = LIBRARY, coords: list = COORDS, primers: list = PRIMERS,
        fill: bool = FILL, parallel: str = PARALLEL):
    """
    Run the vectoring step.
    Generate a vector encoding mutations for each read
    (or read pair, if paired-end).

    FASTA (path): reference sequence(s) in FASTA format

    BAM_FILES (paths): list of one or more alignment files (space separated)
                       in BAM format
    """

    # read library
    if library:
        add_coords_from_library(library, coords, fasta, out_dir)
                        
    # Compute mutation vectors for each BAM file
    bam_files = [os.path.join(bam_dir, bam_file)
                 for bam_dir in bam_dirs
                 for bam_file in os.listdir(bam_dir)
                 if bam_file.endswith(BAM_EXT)]
    writers = VectorWriterSpawner(out_dir, fasta, bam_files, coords,
                                  encode_primers(primers), fill, parallel)
    writers.profile()
