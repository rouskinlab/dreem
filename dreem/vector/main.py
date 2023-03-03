import logging
import os
from typing import Iterable

from click import command, pass_obj
import pandas as pd

from ..align.reads import BamIndexer
from ..util.cli import (DreemCommandName, dreem_command,
                        opt_fasta, opt_bamf, opt_bamd,
                        opt_library, opt_cfill,
                        opt_coords, opt_primers, opt_primer_gap,
                        opt_out_dir, opt_temp_dir,
                        opt_phred_enc, opt_min_phred,
                        opt_strict_pairs, opt_ambindel, opt_batch_size,
                        opt_parallel, opt_max_procs,
                        opt_rerun, opt_resume, opt_save_temp)
from ..util.docdef import autodef, autodoc
from ..util.path import BAM_EXT, OneRefAlignmentInFilePath, RefsetSeqInFilePath
from ..util.seq import DNA
from ..vector.profile import generate_profiles, get_writers


def add_coords_from_library(library_path: str,
                            coords: list[tuple[str, int, int]]):
    try:
        library = pd.read_csv(library_path)
        for ref, end5, end3 in zip(library["construct"],
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
                if coord in coords:
                    logging.warning(f"Skipping duplicate coordinates: {coord}")
                else:
                    coords.append(coord)
    except (FileNotFoundError, KeyError, ValueError) as error:
        logging.error(f"Failed to add coordinates from {library_path} "
                      f"with the following error: {error}")


def encode_primers(primers: Iterable[tuple[str, str, str]]):
    for ref, fwd, rev in primers:
        try:
            yield ref, DNA(fwd.encode()), DNA(rev.encode())
        except ValueError as error:
            logging.error(f"Failed to add primer pair {ref, fwd, rev} "
                          f"with the following error: {error}")


def list_bam_paths(bamf: tuple[str, ...], bamd: tuple[str, ...]):
    bam_paths: dict[str, OneRefAlignmentInFilePath] = dict()

    def add_bam_path(bam):
        bam_path = str(bam)
        if not os.path.isfile(bam_path):
            logging.error(f"Skipping non-existant BAM file: {bam_path}")
            return
        if not bam_path.endswith(BAM_EXT):
            logging.error(f"Skipping non-BAM-formatted file: {bam_path}")
            return
        if bam_path in bam_paths:
            logging.warning(f"Skipping duplicate BAM file: {bam_path}")
            return
        bam_paths[bam_path] = OneRefAlignmentInFilePath.parse(bam_path)

    for bam_file in bamf:
        add_bam_path(bam_file)
    for bam_dir in bamd:
        if not os.path.isdir(bam_dir):
            logging.error(f"No such BAM directory: {bam_dir}")
            continue
        for bam_file in os.listdir(bam_dir):
            if bam_file.endswith(BAM_EXT):
                add_bam_path(os.path.join(bam_dir, bam_file))

    return list(bam_paths.values())


@command(DreemCommandName.VECTOR.value, params=[
    # Input files,
    opt_fasta,
    opt_bamf,
    opt_bamd,
    # Regions
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
    opt_ambindel,
    opt_strict_pairs,
    opt_batch_size,
    # Parallelization
    opt_parallel,
    opt_max_procs,
    # File generation
    opt_rerun,
    opt_resume,
    opt_save_temp,
])
# Pass context object
@pass_obj
# Turn into DREEM command
@dreem_command(imports=("fasta", "bamf"),
               exports="mp_report")
def cli(*args, **kwargs):
    """ Hook the command line interface to the ```run``` function. """
    return run(*args, **kwargs)


@autodoc()
@autodef()
def run(fasta: str,
        /, *,
        bamf: tuple[str],
        bamd: tuple[str],
        coords: tuple[tuple[str, int, int], ...],
        primers: tuple[tuple[str, str, str], ...],
        primer_gap: int,
        cfill: bool,
        library: str,
        batch_size: float,
        max_procs: int,
        **kwargs):
    """ Run the vectoring step. Generate a vector encoding mutations for
    each read (or read pair, if paired-end). """

    # Convert reference file string to path.
    refset_path = RefsetSeqInFilePath.parse(fasta)
    # List out all the paths to BAM files.
    bam_paths = list_bam_paths(bamf, bamd)
    # Convert given primer sequences (str) to DNA objects.
    primers = tuple(encode_primers(primers))
    # Convert batch_size from mebibytes (2^20 = 1048576 bytes) to bytes.
    bytes_per_batch = round(batch_size * 1048576)

    # Index every BAM file.
    for bam_path in bam_paths:
        BamIndexer(xam=bam_path, num_cpus=max_procs, resume=True).run()

    # If a library file is given, add coordinates from the file.
    if library:
        coords = list(coords)
        add_coords_from_library(library_path=library,
                                coords=coords)
        coords = tuple(coords)

    # Compute mutation mut_vectors for each BAM file.
    writers = get_writers(refset_path, bam_paths,
                          coords=coords,
                          primers=primers,
                          primer_gap=primer_gap,
                          cfill=cfill)
    return generate_profiles(writers,
                             batch_size=bytes_per_batch,
                             max_procs=max_procs,
                             **kwargs)
