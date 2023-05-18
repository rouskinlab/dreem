from concurrent.futures import Future, ProcessPoolExecutor
from functools import partial
from logging import getLogger
from pathlib import Path

from click import command
from .krun import cluster_sect
from ..util import docdef
from ..util.cli import (opt_out_dir, opt_mv_file, opt_parallel,
                        opt_max_procs, opt_coords, opt_primers, opt_primer_gap,
                        opt_library, opt_info_thresh,
                        opt_max_clusters, opt_num_runs, opt_signal_thresh,
                        opt_exclude_gu, opt_include_del, opt_exclude_polya,
                        opt_min_iter, opt_max_iter, opt_convergence_cutoff,
                        opt_min_reads, opt_include_ins, opt_min_gap,
                        opt_max_muts_per_read)
from ..util.parallel import get_num_parallel
from ..util.sect import encode_primers
from ..vector.load import open_sections


logger = getLogger(__name__)


params = [
    # Input/output paths
    opt_mv_file,
    opt_out_dir,
    # Sections
    opt_coords,
    opt_primers,
    opt_primer_gap,
    opt_library,
    # Parallelization
    opt_max_procs,
    opt_parallel,
    # Clustering options
    opt_max_clusters,
    opt_num_runs,
    opt_info_thresh,
    opt_signal_thresh,
    opt_exclude_polya,
    opt_exclude_gu,
    opt_include_del,
    opt_include_ins,
    opt_min_gap,
    opt_max_muts_per_read,
    opt_min_iter,
    opt_max_iter,
    opt_convergence_cutoff,
    opt_min_reads
]


@command("cluster", params=params)
def cli(*args, **kwargs):
    return run(*args, **kwargs)


@docdef.auto()
def run(mv_file: tuple[str, ...], *,
        out_dir: str,
        # Sections
        library: str,
        coords: tuple[tuple[str, int, int], ...],
        primers: tuple[tuple[str, str, str], ...],
        primer_gap: int,
        max_procs: int,
        parallel: bool,
        # Clustering options
        max_clusters: int,
        num_runs: int,
        signal_thresh: float,
        exclude_gu: bool,
        include_del: bool,
        include_ins: bool,
        exclude_polya: int,
        max_muts_per_read: int,
        min_gap: int,
        info_thresh: float,
        min_iter: int,
        max_iter: int,
        convergence_cutoff: float,
        min_reads: int) -> list[Path]:
    """ Run the clustering module. """
    # Open all vector reports and get the sections for each.
    reports, sections = open_sections(map(Path, mv_file),
                                      coords=coords,
                                      primers=encode_primers(primers),
                                      primer_gap=primer_gap,
                                      library=(Path(library) if library
                                               else None))
    # Determine how to parallelize clustering.
    n_tasks, n_procs_per_task = get_num_parallel(sections.count, max_procs,
                                                 parallel=parallel,
                                                 hybrid=False)
    # Run EM clustering on every section of every set of bit vectors.
    cluster_func = partial(cluster_sect,
                           exclude_gu=exclude_gu, include_del=include_del,
                           include_ins=include_ins, exclude_polya=exclude_polya,
                           max_muts_per_read=max_muts_per_read,
                           min_gap=min_gap, min_reads=min_reads,
                           info_thresh=info_thresh, signal_thresh=signal_thresh,
                           max_clusters=max_clusters, n_runs=num_runs,
                           conv_thresh=convergence_cutoff,
                           min_iter=min_iter, max_iter=max_iter,
                           max_procs=n_procs_per_task, out_dir=Path(out_dir))
    if n_tasks > 1:
        # Run multiple sections of bit vectors in parallel.
        with ProcessPoolExecutor(max_workers=n_tasks) as pool:
            futures: list[Future] = list()
            for report in reports:
                for sect in sections.list(report.ref):
                    futures.append(pool.submit(cluster_func, report, sect))
                    logger.debug(f"Submitted EM clustering for {report} {sect}")
            cluster_files: list[Path] = [future.result() for future in futures]
    else:
        # Run each section of bit vectors one at a time.
        cluster_files: list[Path] = list()
        for report in reports:
            for sect in sections.list(report.ref):
                cluster_files.append(cluster_func(report, sect))
    # Remove any None values (indicating that clustering failed).
    return list(filter(None, cluster_files))
