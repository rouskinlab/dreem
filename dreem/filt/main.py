from concurrent.futures import Future, ProcessPoolExecutor
from functools import partial
from logging import getLogger
from pathlib import Path

from click import command

from .filt import filter_sect
from ..mvec.load import open_sections
from ..util import docdef
from ..util.bit import BitCaller
from ..util.cli import (opt_mvec,
                        opt_coords, opt_primers, opt_primer_gap, opt_library,
                        opt_count_del, opt_count_ins, opt_discount_mut,
                        opt_exclude_polya, opt_exclude_gu, opt_exclude_pos,
                        opt_min_finfo_read, opt_max_fmut_read, opt_min_mut_gap,
                        opt_min_ninfo_pos, opt_max_fmut_pos,
                        opt_max_procs, opt_parallel)
from ..util.parallel import get_num_parallel
from ..util.sect import encode_primers

logger = getLogger(__name__)

params = [
    # Input/output paths
    opt_mvec,
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


@command("filter", params=params)
def cli(*args, **kwargs):
    return run(*args, **kwargs)


@docdef.auto()
def run(mv_file: tuple[str, ...], *,
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
    # Open all vector reports and get the sections for each.
    loaders, sections = open_sections(map(Path, mv_file),
                                      coords=coords,
                                      primers=encode_primers(primers),
                                      primer_gap=primer_gap,
                                      library=(Path(library) if library
                                               else None))
    # Determine how to parallelize filtering.
    n_tasks, n_procs_per_task = get_num_parallel(sections.count, max_procs,
                                                 parallel=parallel,
                                                 hybrid=False)
    bit_caller = BitCaller.from_counts(count_del=count_del,
                                       count_ins=count_ins,
                                       discount=discount_mut)
    # Filter every section of every set of bit vectors.
    filter_func = partial(filter_sect,
                          bit_caller=bit_caller,
                          exclude_polya=exclude_polya,
                          exclude_gu=exclude_gu,
                          exclude_pos=exclude_pos,
                          min_finfo_read=min_finfo_read,
                          max_fmut_read=max_fmut_read,
                          min_mut_gap=min_mut_gap,
                          min_ninfo_pos=min_ninfo_pos,
                          max_fmut_pos=max_fmut_pos)
    if n_tasks > 1:
        # Run multiple sections of bit vectors in parallel.
        with ProcessPoolExecutor(max_workers=n_tasks) as pool:
            futures: list[Future] = list()
            for loader in loaders:
                for section in sections.list(loader.ref):
                    futures.append(pool.submit(filter_func, loader, section))
                    logger.debug(f"Submitted filtering for {loader} {section}")
            filter_files: list[Path] = [future.result() for future in futures]
    else:
        # Filter each section of bit vectors one at a time.
        filter_files: list[Path] = list()
        for loader in loaders:
            for section in sections.list(loader.ref):
                filter_files.append(filter_func(loader, section))
    # Remove any None values (indicating that clustering failed).
    return list(filter(None, filter_files))
