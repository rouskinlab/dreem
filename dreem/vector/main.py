from logging import getLogger
import pathlib

from click import command

from ..util import docdef, path
from ..util.cli import (opt_fasta, opt_bamf, opt_bamd,
                        opt_out_dir, opt_temp_dir,
                        opt_phred_enc, opt_min_phred,
                        opt_ambid, opt_batch_size,
                        opt_parallel, opt_max_procs,
                        opt_rerun, opt_save_temp)
from ..util.parallel import lock_temp_dir
from ..vector.write import write, get_writers


logger = getLogger(__name__)


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
    # Output directories
    opt_out_dir,
    opt_temp_dir,
    # SAM options
    opt_phred_enc,
    opt_min_phred,
    # Vectoring options
    opt_ambid,
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


@lock_temp_dir
@docdef.auto()
def run(fasta: str,
        *,
        bamf: tuple[str],
        bamd: tuple[str],
        out_dir: str,
        temp_dir: str,
        phred_enc: int,
        min_phred: int,
        ambid: bool,
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

    # Create an object to write mutation vectors for each BAM file.
    writers = get_writers(pathlib.Path(fasta), list_bam_paths(bamf, bamd))

    # Compute and write mutation vectors for each BAM file.
    profiles = write(writers=writers,
                     out_dir=pathlib.Path(out_dir),
                     temp_dir=pathlib.Path(temp_dir),
                     phred_enc=phred_enc,
                     min_phred=min_phred,
                     ambid=ambid,
                     batch_size=batch_size,
                     max_procs=max_procs,
                     parallel=parallel,
                     rerun=rerun,
                     save_temp=save_temp)
    return profiles
