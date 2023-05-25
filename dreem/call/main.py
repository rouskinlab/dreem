from logging import getLogger
from pathlib import Path

from click import command

from .filt import filter_sect
from ..relate.load import open_sections
from ..core import docdef, path
from dreem.core.bit import BitCaller
from ..core.cli import (opt_rel,
                        opt_coords, opt_primers, opt_primer_gap, opt_library,
                        opt_count_del, opt_count_ins, opt_discount_mut,
                        opt_exclude_polya, opt_exclude_gu, opt_exclude_pos,
                        opt_min_finfo_read, opt_max_fmut_read, opt_min_mut_gap,
                        opt_min_ninfo_pos, opt_max_fmut_pos,
                        opt_max_procs, opt_parallel)
from ..core.parallel import dispatch
from ..core.sect import encode_primers

logger = getLogger(__name__)

params = [
    # Input/output paths
    opt_rel,
    # Sections
    opt_coords, opt_primers, opt_primer_gap, opt_library,
    # Mutation counting
    opt_count_del, opt_count_ins, opt_discount_mut,
    # Filtering
    opt_exclude_polya, opt_exclude_gu, opt_exclude_pos,
    opt_min_finfo_read, opt_max_fmut_read, opt_min_mut_gap,
    opt_min_ninfo_pos, opt_max_fmut_pos,
    # Parallelization
    opt_max_procs, opt_parallel,
]


@command(path.MOD_CALL, params=params)
def cli(*args, **kwargs):
    return run(*args, **kwargs)


@docdef.auto()
def run(rel: tuple[str, ...], *,
        # Sections
        coords: tuple[tuple[str, int, int], ...],
        primers: tuple[tuple[str, str, str], ...],
        primer_gap: int,
        library: str,
        # Mutation counting
        count_del: bool,
        count_ins: bool,
        discount_mut: tuple[str, ...],
        # Filtering
        exclude_polya: int,
        exclude_gu: bool,
        exclude_pos: tuple[tuple[str, int], ...],
        min_finfo_read: float,
        max_fmut_read: int,
        min_mut_gap: int,
        min_ninfo_pos: int,
        max_fmut_pos: float,
        # Parallelization
        max_procs: int,
        parallel: bool) -> list[Path]:
    """ Run the filtering module. """
    # Open all relation vector loaders and get the sections for each.
    loaders, sections = open_sections(map(Path, rel),
                                      coords=coords,
                                      primers=encode_primers(primers),
                                      primer_gap=primer_gap,
                                      library=(Path(library) if library
                                               else None))
    # List the relation loaders and their sections.
    args = [(loader, section) for loader in loaders
            for section in sections.list(loader.ref)]
    # Define the keyword arguments.
    kwargs = dict(bit_caller=BitCaller.from_counts(count_del=count_del,
                                                   count_ins=count_ins,
                                                   discount=discount_mut),
                  exclude_polya=exclude_polya,
                  exclude_gu=exclude_gu,
                  exclude_pos=exclude_pos,
                  min_finfo_read=min_finfo_read,
                  max_fmut_read=max_fmut_read,
                  min_mut_gap=min_mut_gap,
                  min_ninfo_pos=min_ninfo_pos,
                  max_fmut_pos=max_fmut_pos)
    # Call the mutations and filter the mutation vectors.
    reports = dispatch(filter_sect, max_procs=max_procs, parallel=parallel,
                       pass_n_procs=False, args=args, kwargs=kwargs)
    return reports
