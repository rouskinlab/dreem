import os
from typing import List, Tuple

import click
import pandas as pd

from dreem.util.util import DNA
from dreem.util.cli import opto_out_dir, opti_library, opti_coords, opti_primers, opti_fill, opti_parallel, argi_fasta, argi_bams
from dreem.vector.mprofile import VectorWriterSpawner
from dreem.util.files_sanity import check_library

def add_coords_from_library(library_path: str,
                            coords: List[Tuple[str, int, int]],
                            fasta: str):
    library = check_library(pd.read_csv(library_path), fasta)
    for row in library.index:
        construct = os.path.splitext(library.loc[row, "construct"])[0]
        # Convert from inclusive 0-indexed "section_start"
        # to inclusive 1-indexed "first"
        first = int(library.loc[row, "section_start"]) + 1
        # Convert from exclusive 0-indexed "section_end"
        # to inclusive 1-indexed "last"
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


@click.command()
@opti_library
@opti_coords
@opti_primers
@opti_fill
@opti_parallel
@opto_out_dir
@argi_fasta
@argi_bams
def run(out_dir: str, fasta: str, bam_files: List[str],
        library: str, coords: list, primers: list,
        fill: bool, parallel: str):
    """
    Run the vectoring step.
    Generate a vector encoding mutations for each read
    (or read pair, if paired-end).

    FASTA (path): reference sequence(s) in FASTA format

    BAM_FILES (paths): list of one or more alignment files (space separated)
                       in BAM format
    """

    # Create the folders
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    
    # read library
    if library:
        add_coords_from_library(library, coords, fasta)
    
    # encode coordinates and primers
    coords, primers = encode_coords(coords, primers)
    
    # Compute mutation vectors for each BAM file
    writers = VectorWriterSpawner(out_dir, fasta, bam_files, coords, primers,
                                  fill, parallel)
    writers.gen_mut_profiles()


if __name__ == "__main__":
    run()
