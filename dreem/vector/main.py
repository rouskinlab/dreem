import os
from typing import List, Tuple

import pandas as pd
import subprocess
import warnings

from dreem.util.util import DNA, run_cmd
from dreem.util.cli import OUT_DIR, LIBRARY, COORDS, PRIMERS, FILL, PARALLEL
from dreem.util.path import MOD_VEC
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


def encode_coords(coords: List[Tuple[str, int, int]],
                  primers: List[Tuple[str, str, str]]):
    coords_bytes = [(ref.encode(), first, last)
                    for ref, first, last in coords]
    primers_bytes = [(ref.encode(), DNA(fwd.encode()), DNA(rev.encode()))
                     for ref, fwd, rev in primers]
    return coords_bytes, primers_bytes


def run(fasta: str, bam_files: List[str], out_dir: str = OUT_DIR,
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

    # Add "vectoring" to outdir
    out_dir = os.path.join(out_dir, MOD_VEC)

    # Create the directory
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    
    # read library
    if library:
        add_coords_from_library(library, coords, fasta, out_dir)
    
    # encode coordinates and primers
    coords, primers = encode_coords(coords, primers)
                        
    # Compute mutation vectors for each BAM file
    writers = VectorWriterSpawner(out_dir, fasta, bam_files, coords, primers,
                                  fill, parallel)
    writers.gen_mut_profiles()
