import os

import pandas as pd

from dreem.util.seq import DNA
from dreem.util.path import BAM_EXT
from dreem.vector.mprofile import VectorWriterSpawner
from dreem.util.files_sanity import check_library


def add_coords_from_library(library_path: str,
                            coords: list[tuple[str, int, int]],
                            fasta: str,
                            out_dir: str):
    library = check_library(pd.read_csv(library_path), fasta, out_dir)
    for construct, end5, end3 in zip(library["construct"],
                                     library["section_start"],
                                     library["section_end"], strict=True):
        coord = (str(os.path.splitext(construct)[0]), int(end5), int(end3))
        if coord not in coords:
            coords.append(coord)


def encode_primers(primers: list[tuple[str, str, str]]):
    return [(ref, DNA(fwd.encode()), DNA(rev.encode()))
            for ref, fwd, rev in primers]


def run(out_dir: str,
        temp_dir: str,
        bams: tuple[str],
        fasta: str,
        coords: tuple[tuple[str, int, int], ...],
        primers: tuple[tuple[str, str, str], ...],
        primer_gap: int,
        spanall: bool,
        library: str,
        parallel: bool,
        max_procs: int,
        phred_enc: int,
        min_phred: int,
        rerun: bool,
        save_temp: bool,
        resume: bool):
    """
    Run the vectoring step.
    Generate a vector encoding mutations for each read
    (or read pair, if paired-end).

    FASTA (path): reference sequence(s) in FASTA format

    BAM_FILES (paths): list of one or more alignment files (space separated)
                       in BAM format
    """

    # Change data types
    bams = list(bams)
    coords = list(coords)
    primers = encode_primers(list(primers))

    # If a library file is given, add coordinates from the file.
    if library:
        add_coords_from_library(library_path=library,
                                coords=coords,
                                fasta=fasta,
                                out_dir=out_dir)
                        
    # Compute mutation vectors for each BAM file
    writers = VectorWriterSpawner(out_dir=out_dir,
                                  temp_dir=temp_dir,
                                  fasta=fasta,
                                  bam_files=bams,
                                  coords=coords,
                                  primers=primers,
                                  primer_gap=primer_gap,
                                  spanall=spanall,
                                  parallel=parallel,
                                  max_procs=max_procs,
                                  min_phred=min_phred,
                                  phred_enc=phred_enc,
                                  rerun=rerun)
    return writers.generate_profiles(save_temp, resume)
