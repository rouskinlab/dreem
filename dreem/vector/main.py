from logging import getLogger
import pathlib
from typing import Iterable

from click import command
import pandas as pd

from ..util import docdef, path
from ..util.cli import (opt_fasta, opt_bamf, opt_bamd,
                        opt_library, opt_autosect,
                        opt_coords, opt_primers, opt_primer_gap,
                        opt_out_dir, opt_temp_dir,
                        opt_phred_enc, opt_min_phred,
                        opt_strict_pairs, opt_ambid, opt_batch_size,
                        opt_parallel, opt_max_procs,
                        opt_rerun, opt_save_temp)
from ..util.parallel import lock_output
from ..util.seq import DNA
from ..vector.profile import generate_profiles, get_writers


logger = getLogger(__name__)


def add_coords_from_library(library_path: str,
                            initial_coords: tuple[tuple[str, int, int]]):
    library_coords = list()
    try:
        library = pd.read_csv(library_path)
        for ref, end5, end3 in zip(library["reference"],
                                   library["section_start"],
                                   library["section_end"], strict=True):
            try:
                coord = (str(pathlib.Path(ref).with_suffix("")),
                         int(end5),
                         int(end3))
            except ValueError as error:
                logger.error(f"Failed to add coordinates {ref, end5, end3} "
                              f"with the following error: {error}")
            else:
                if coord in initial_coords or coord in library_coords:
                    logger.warning(f"Skipping duplicate coordinates: {coord}")
                else:
                    library_coords.append(coord)
    except (FileNotFoundError, KeyError, ValueError) as error:
        logger.critical(f"Failed to add coordinates from {library_path} "
                      f"with the following error: {error}")
    return initial_coords + tuple(library_coords)


def encode_primers(primers: Iterable[tuple[str, str, str]]):
    for ref, fwd, rev in primers:
        try:
            yield ref, DNA(fwd.encode()), DNA(rev.encode())
        except ValueError as error:
            logger.error(f"Failed to add primer pair {ref, fwd, rev} "
                          f"with the following error: {error}")


def list_bam_paths(bamf: tuple[str, ...], bamd: tuple[str, ...]):
    bam_files: set[pathlib.Path] = set()

    def add_bam_file(file: pathlib.Path):
        if not file.is_file():
            logger.critical(f"Skipping non-existant BAM file: {file}")
            return
        if file.suffix != path.BAM_EXT:
            logger.critical(f"Skipping non-BAM-formatted file: {file}")
            return
        if file in bam_files:
            logger.warning(f"Skipping duplicate BAM file: {file}")
            return
        bam_files.add(file)

    for bam_file in bamf:
        add_bam_file(pathlib.Path(bam_file))
    for bam_dir in bamd:
        dpath = pathlib.Path(bam_dir)
        if not dpath.is_dir():
            logger.critical(f"Skipping non-existant BAM directory: {dpath}")
            continue
        for bamd_file in dpath.iterdir():
            add_bam_file(bamd_file)

    return list(bam_files)


# Parameters for command line interface
params = [
    # Input files
    opt_fasta,
    opt_bamf,
    opt_bamd,
    # Sections
    opt_library,
    opt_autosect,
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
    opt_save_temp,
]


@command("vector", params=params)
def cli(**kwargs):
    """ Hook the command line interface to the ```run``` function. """
    return run(**kwargs)


@lock_output
@docdef.auto()
def run(fasta: str,
        *,
        bamf: tuple[str],
        bamd: tuple[str],
        library: str,
        autosect: bool,
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
        save_temp: bool):
    """ Run the vectoring step. Generate a vector encoding mutations for
    each read (or read pair, if paired-end). """

    if not fasta:
        logger.critical("No FASTA file was given to vectoring")
        return ()

    # If a library file is given, add coordinates from the file.
    if library:
        coords_lib = add_coords_from_library(library_path=library,
                                             initial_coords=coords)
    else:
        coords_lib = coords

    # List all BAM files among the given files and directories.
    bams = list_bam_paths(bamf, bamd)

    # Create an object to write mutation vectors for each BAM file.
    writers, temp_bams = get_writers(pathlib.Path(fasta), bams,
                                     temp_dir=pathlib.Path(temp_dir),
                                     n_procs=max_procs,
                                     coords=coords_lib,
                                     primers=tuple(encode_primers(primers)),
                                     primer_gap=primer_gap,
                                     autosect=autosect)

    try:
        # Compute and write mutation vectors for each BAM file.
        profiles = generate_profiles(writers=writers,
                                     out_dir=pathlib.Path(out_dir),
                                     temp_dir=pathlib.Path(temp_dir),
                                     phred_enc=phred_enc,
                                     min_phred=min_phred,
                                     ambid=ambid,
                                     strict_pairs=strict_pairs,
                                     batch_size=batch_size,
                                     max_procs=max_procs,
                                     parallel=parallel,
                                     rerun=rerun,
                                     save_temp=save_temp)
        return profiles
    finally:
        if not save_temp:
            # Always delete the temporary BAM and index files, whether
            # vectoring completed normally or raised an error.
            for temp_bam in temp_bams:
                temp_bam.unlink(missing_ok=True)
