import os
from typing import List, Tuple

import pandas as pd

from dreem.util.util import DNA
from dreem.util.path import BAM_EXT
from dreem.vector.mprofile import VectorWriterSpawner
from dreem.util.files_sanity import check_library


def add_coords_from_library(library_path: str,
                            coords: List[Tuple[str, int, int]],
                            fasta: str,
                            out_dir: str):
    library = check_library(pd.read_csv(library_path), fasta, out_dir)
    for construct, end5, end3 in zip(library["construct"],
                                     library["section_start"],
                                     library["section_end"], strict=True):
        coord = (str(os.path.splitext(construct)[0]), int(end5), int(end3))
        if coord not in coords:
            coords.append(coord)


def encode_primers(primers: List[Tuple[str, str, str]]):
    return [(ref, DNA(fwd.encode()), DNA(rev.encode()))
            for ref, fwd, rev in primers]


def run(fasta: str, bam_dirs: list[str], phred_enc: int, min_phred: int,
        library: str, coords: list, primers: list, fill: bool, parallel: str,
        max_cpus: int, top_dir: str, rerun: bool):
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
        add_coords_from_library(library_path=library,
                                coords=coords,
                                fasta=fasta,
                                out_dir=top_dir)
                        
    # Compute mutation vectors for each BAM file
    bam_files = [os.path.join(bam_dir, bam_file)
                 for bam_dir in set(bam_dirs)
                 for bam_file in os.listdir(bam_dir)
                 if bam_file.endswith(BAM_EXT)]
    writers = VectorWriterSpawner(top_dir=top_dir,
                                  fasta=fasta,
                                  bam_files=bam_files,
                                  coords=coords,
                                  primers=encode_primers(primers),
                                  fill=fill,
                                  parallel=parallel,
                                  max_cpus=max_cpus,
                                  min_phred=min_phred,
                                  phred_enc=phred_enc,
                                  rerun=rerun)
    return writers.generate_profiles()
