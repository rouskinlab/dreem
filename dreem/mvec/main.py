from logging import getLogger
from pathlib import Path

from click import command

from ..util import docdef, path
from ..util.cli import (opt_fasta, opt_bam,
                        opt_out_dir, opt_temp_dir,
                        opt_phred_enc, opt_min_phred,
                        opt_ambid, opt_batch_size,
                        opt_parallel, opt_max_procs,
                        opt_rerun, opt_save_temp)
from ..util.parallel import lock_temp_dir
from ..mvec.write import write, get_writers


logger = getLogger(__name__)


# Parameters for command line interface
params = [
    # Input files
    opt_fasta,
    opt_bam,
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
        bam: tuple[str, ...],
        *,
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
        return list()

    # Create an object to write mutation vectors for each BAM file.
    writers = get_writers(Path(fasta),
                          path.find_files_multi(map(Path, bam),
                                                [path.SampSeg, path.XamSeg]))

    # Compute and write mutation vectors for each BAM file.
    profiles = write(writers=writers,
                     out_dir=Path(out_dir),
                     temp_dir=Path(temp_dir),
                     phred_enc=phred_enc,
                     min_phred=min_phred,
                     ambid=ambid,
                     batch_size=batch_size,
                     max_procs=max_procs,
                     parallel=parallel,
                     rerun=rerun,
                     save_temp=save_temp)
    return profiles
