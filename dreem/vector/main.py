import logging
import os
from typing import Iterable

from click import command, pass_obj
import pandas as pd

from ..align.reads import index_bam_file
from ..util.cli import (DreemCommandName, dreem_command,
                        opt_fasta, opt_bamf, opt_bamd,
                        opt_library, opt_cfill,
                        opt_coords, opt_primers, opt_primer_gap,
                        opt_out_dir, opt_temp_dir,
                        opt_phred_enc, opt_min_phred,
                        opt_strict_pairs, opt_ambid, opt_batch_size,
                        opt_parallel, opt_max_procs,
                        opt_rerun, opt_resume, opt_save_temp)
from ..util import docdef
from ..util.path import BAM_EXT
from ..util.seq import DNA
from ..vector.profile import generate_profiles, get_writers


def add_coords_from_library(library_path: str,
                            initial_coords: tuple[tuple[str, int, int]]):
    library_coords = list()
    try:
        library = pd.read_csv(library_path)
        for ref, end5, end3 in zip(library["reference"],
                                   library["section_start"],
                                   library["section_end"], strict=True):
            try:
                coord = (str(os.path.splitext(ref)[0]),
                         int(end5),
                         int(end3))
            except ValueError as error:
                logging.error(f"Failed to add coordinates {ref, end5, end3} "
                              f"with the following error: {error}")
            else:
                if coord in initial_coords or coord in library_coords:
                    logging.warning(f"Skipping duplicate coordinates: {coord}")
                else:
                    library_coords.append(coord)
    except (FileNotFoundError, KeyError, ValueError) as error:
        logging.error(f"Failed to add coordinates from {library_path} "
                      f"with the following error: {error}")
    return initial_coords + tuple(library_coords)


def encode_primers(primers: Iterable[tuple[str, str, str]]):
    for ref, fwd, rev in primers:
        try:
            yield ref, DNA(fwd.encode()), DNA(rev.encode())
        except ValueError as error:
            logging.error(f"Failed to add primer pair {ref, fwd, rev} "
                          f"with the following error: {error}")


def list_bam_paths(bamf: tuple[str, ...], bamd: tuple[str, ...]):
    bam_files: set[str] = set()

    def add_bam_file(file):
        if not os.path.isfile(file):
            logging.error(f"Skipping non-existant BAM file: {file}")
            return
        if os.path.splitext(file)[1] != BAM_EXT:
            logging.error(f"Skipping non-BAM-formatted file: {file}")
            return
        if file in bam_files:
            logging.warning(f"Skipping duplicate BAM file: {file}")
            return
        bam_files.add(file)

    for bam_file in bamf:
        add_bam_file(bam_file)
    for bam_dir in bamd:
        if not os.path.isdir(bam_dir):
            logging.error(f"Skipping non-existant BAM directory: {bam_dir}")
            continue
        for bam_file in os.listdir(bam_dir):
            if os.path.splitext(bam_file)[1] == BAM_EXT:
                add_bam_file(os.path.join(bam_dir, bam_file))

    return list(bam_files)


# Parameters for command line interface
params = [
    # Input files
    opt_fasta,
    opt_bamf,
    opt_bamd,
    # Sections
    opt_library,
    opt_cfill,
    opt_coords,
    opt_primers,
    opt_primer_gap,
    # Output directories
    opt_out_dir,
    opt_temp_dir,
    # SAM options
    opt_phred_enc,
    opt_min_phred,
    # Vectoring options
    opt_ambid,
    opt_strict_pairs,
    opt_batch_size,
    # Parallelization
    opt_max_procs,
    opt_parallel,
    # File generation
    opt_rerun,
    opt_resume,
    opt_save_temp,
]


@command(DreemCommandName.VECTOR.value, params=params)
# Pass context object.
@pass_obj
# Turn into DREEM command.
@dreem_command(imports=("fasta", "bamf"),
               result_key="report")
def cli(**kwargs):
    """ Hook the command line interface to the ```run``` function. """
    return run(**kwargs)


@docdef.auto()
def run(fasta: str,
        *,
        bamf: tuple[str],
        bamd: tuple[str],
        library: str,
        cfill: bool,
        coords: tuple[tuple[str, int, int], ...],
        primers: tuple[tuple[str, str, str], ...],
        primer_gap: int,
        out_dir: str,
        temp_dir: str,
        phred_enc: int,
        min_phred: int,
        ambid: bool,
        strict_pairs: bool,
        batch_size: float,
        max_procs: int,
        parallel: bool,
        rerun: bool,
        resume: bool,
        save_temp: bool):
    """ Run the vectoring step. Generate a vector encoding mutations for
    each read (or read pair, if paired-end). """

    # Index every BAM file.
    bams = list_bam_paths(bamf, bamd)
    for bam in bams:
        index_bam_file(bam)

    # If a library file is given, add coordinates from the file.
    if library:
        coords_lib = add_coords_from_library(library_path=library,
                                             initial_coords=coords)
    else:
        coords_lib = coords

    # Convert the primer sequences from str to DNA.
    primers_enc = tuple(encode_primers(primers))

    # Create an object to write mutation vectors for each BAM file.
    writers = get_writers(fasta,
                          bams,
                          coords=coords_lib,
                          primers=primers_enc,
                          primer_gap=primer_gap,
                          cfill=cfill)

    # Compute and write mutation vectors for each BAM file.
    return generate_profiles(writers=writers,
                             out_dir=out_dir,
                             temp_dir=temp_dir,
                             phred_enc=phred_enc,
                             min_phred=min_phred,
                             ambid=ambid,
                             strict_pairs=strict_pairs,
                             batch_size=batch_size,
                             max_procs=max_procs,
                             parallel=parallel,
                             rerun=rerun,
                             resume=resume,
                             save_temp=save_temp)
