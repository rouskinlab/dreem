import os
from typing import List, Tuple

import pandas as pd

from dreem.util.util import DNA
from dreem.util.cli import FASTA, INPUT_DIR, OUT_DIR, LIBRARY, PARALLEL, COORDS, PRIMERS, FILL
from dreem.vector.mprofile import VectorWriterSpawner
from dreem.util.files_sanity import check_library

def add_coords_from_library(library_path: str,
                            coords: List[Tuple[str, int, int]]):
    library = check_library(pd.read_csv(library_path))
    for row in library.index:
        construct = os.path.splitext(library.loc[row, "construct"])[0]
        first = library.loc[row, "section_start"]
        last = library.loc[row, "section_end"]
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


def run(out_dir: str, fasta: str, bam_files: List[str],
        library:str=LIBRARY, coords:list=COORDS, primers:list=PRIMERS,
        fill:bool=FILL, parallel:str=PARALLEL):
    """Run the vectoring pipeline.

    Turns each bam file into a vector file and outputs them in the directory `out_dir`.

    Parameters from args:
    -----------------------
    out_dir: str
        Path to the output folder.
    fasta: str
        Path to the reference FASTA file.
    bam_files: List[str]
        List of path to each BAM file to process.
    coords: tuple
        coordinates for reference: '-c ref-name first last'
    primers: tuple
        primers for reference: '-p ref-name fwd rev'
    fill: bool
        Fill in coordinates of reference sequences for which neither coordinates nor primers were given (default: no).
    library: str
        Name of the library. Default is None.
    parallel: str
        Parallelize the processing of mutational PROFILES or READS within each profile, turn parallelization OFF, or AUTO matically choose the parallelization method (default: auto).
    
    Returns
    -------
    1 if successful, 0 otherwise.

    """

    # Create the folders
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    
    # read library
    if library:
        add_coords_from_library(library, coords)
    
    # encode coordinates and primers
    coords, primers = encode_coords(coords, primers)
    
    # Compute mutation vectors for each BAM file
    writers = VectorWriterSpawner(out_dir, fasta, bam_files, coords, primers,
                                  fill, parallel)
    writers.gen_mut_profiles()
