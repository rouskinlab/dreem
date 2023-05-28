from logging import getLogger
from pathlib import Path

from click import command

from .rnastructure import fold
from ..core import docdef, path
from ..core.cli import (opt_temp_dir, opt_save_temp, opt_table,
                        opt_fasta, opt_library,
                        opt_coords, opt_primers, opt_primer_gap,
                        opt_dms_quantile,
                        opt_rnastructure_use_temp,
                        opt_rnastructure_fold_args, opt_rnastructure_use_dms,
                        opt_rnastructure_dms_min_unpaired_value,
                        opt_rnastructure_dms_max_paired_value,
                        opt_rnastructure_deltag_ensemble,
                        opt_rnastructure_probability,
                        opt_max_procs, opt_parallel, opt_rerun)
from ..core.dependencies import check_rnastructure_exists
from ..core.parallel import dispatch, lock_temp_dir
from ..core.rna import RnaProfile
from ..core.sect import RefSections, Section, encode_primers
from ..core.seq import parse_fasta
from ..table.io import load, BitVecPosTableLoader, ClusterPosTableLoader

logger = getLogger(__name__)

params = [
    opt_table,
    opt_fasta,
    opt_library,
    opt_coords,
    opt_primers,
    opt_primer_gap,
    opt_dms_quantile,
    # opt_rnastructure_use_temp,
    # opt_rnastructure_fold_args,
    # opt_rnastructure_use_dms,
    # opt_rnastructure_dms_min_unpaired_value,
    # opt_rnastructure_dms_max_paired_value,
    # opt_rnastructure_deltag_ensemble,
    # opt_rnastructure_probability,
    opt_temp_dir,
    opt_save_temp,
    opt_max_procs,
    opt_parallel,
    opt_rerun,
]


@command(path.MOD_FOLD, params=params)
def cli(**kwargs):
    return run(**kwargs)


@lock_temp_dir
@docdef.auto()
def run(table: tuple[str, ...],
        *,
        fasta: str,
        library: str | None,
        coords: tuple[tuple[str, int, int], ...],
        primers: tuple[tuple[str, str, str], ...],
        primer_gap: int,
        dms_quantile: float,
        temp_dir: str,
        save_temp: bool,
        max_procs: int,
        parallel: bool,
        rerun: bool):
    """
    Run the fold module.
    """
    check_rnastructure_exists()
    # Get the sections for every reference sequence.
    ref_sections = RefSections(parse_fasta(Path(fasta)),
                               library=Path(library) if library else None,
                               coords=coords,
                               primers=encode_primers(primers),
                               primer_gap=primer_gap)
    # Initialize the table loaders.
    loaders = dispatch(load, max_procs, parallel,
                       args=[(Path(file),) for file in table],
                       pass_n_procs=False)
    # Fold the RNA profiles.
    return dispatch(fold_rna, max_procs, parallel,
                    args=[(loader, ref_sections.list(loader.ref))
                          for loader in loaders],
                    kwargs=dict(temp_dir=Path(temp_dir), save_temp=save_temp,
                                dms_quantile=dms_quantile, rerun=rerun),
                    pass_n_procs=True)


def fold_rna(loader: BitVecPosTableLoader | ClusterPosTableLoader,
             sections: list[Section], n_procs: int, **kwargs):
    """ Fold an RNA molecule from one table of reactivities. """
    return dispatch(fold_profile, n_procs, parallel=True,
                    args=[(profile,)
                          for profile in loader.iter_profiles(sections)],
                    kwargs=dict(out_dir=loader.out_dir, **kwargs),
                    pass_n_procs=False)


def fold_profile(rna: RnaProfile, **kwargs):
    return fold(rna, **kwargs)
