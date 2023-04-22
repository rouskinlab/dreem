from pathlib import Path

from click import command
import pandas as pd

from . import tasks
from .bitvector import BitVector
from .em import EMclustering
from ..util import docdef, path
from ..util.cli import (opt_out_dir, opt_mv_file,
                        opt_max_procs, opt_coords, opt_primers, opt_primer_gap,
                        opt_library, opt_info_thresh,
                        opt_max_clusters, opt_num_runs, opt_signal_thresh,
                        opt_include_gu, opt_include_del, opt_max_polya,
                        opt_min_iter, opt_max_iter, opt_convergence_cutoff,
                        opt_min_reads, opt_include_ins, opt_min_mut_dist,
                        opt_max_muts_per_read)
from ..util.sect import RefSections, encode_primers
from ..vector.load import open_reports

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
    # Clustering options
    opt_max_clusters,
    opt_num_runs,
    opt_info_thresh,
    opt_signal_thresh,
    opt_max_polya,
    opt_include_gu,
    opt_include_del,
    opt_include_ins,
    opt_min_mut_dist,
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
        # Clustering options
        max_clusters: int,
        num_runs: int,
        signal_thresh: float,
        include_gu: bool,
        include_del: bool,
        include_ins: bool,
        max_polya: int,
        max_muts_per_read: int,
        min_mut_dist: int,
        info_thresh: float,
        min_iter: int,
        max_iter: int,
        convergence_cutoff: float,
        min_reads: int):
    """ Run the clustering module. """
    out_path = Path(out_dir)
    # Open all vector reports and get the sections for each.
    reports = open_reports(map(Path, mv_file))
    sections = RefSections({(rep.ref, rep.seq) for rep in reports},
                           coords=coords,
                           primers=encode_primers(primers),
                           primer_gap=primer_gap,
                           library=Path(library))

    # Get the bitvector files in the input directory and all of its subdirectories
    cluster_out_files = list()
    for report in reports:
        for sect in sections.list(report.ref):
            # Run the clustering algorithm
            bitvector = BitVector(report, end5=sect.end5, end3=sect.end3,
                                  max_polya=max_polya,
                                  include_gu=include_gu,
                                  include_del=include_del,
                                  include_ins=include_ins,
                                  max_muts_per_read=max_muts_per_read,
                                  min_mut_dist=min_mut_dist,
                                  min_reads=min_reads,
                                  info_thresh=info_thresh,
                                  signal_thresh=signal_thresh)
            bitvector.publish_preprocessing_report(out_path)
            clusters = tasks.run(bitvector, max_clusters, num_runs, min_iter,
                                 max_iter, convergence_cutoff, max_procs)

            # Compute the likelihood of each read for each cluster
            columns = [(k, c + 1) for k in clusters for c in range(k)]
            reads_best_cluster = pd.DataFrame(dtype=float,
                                              index=pd.Index(bitvector.read_names, name="Read Name"),
                                              columns=pd.MultiIndex.from_tuples(columns,
                                                                                names=["K", "i"]))
            for k, clusters_k in clusters.items():
                em = EMclustering(bitvector, k, min_iter, max_iter, convergence_cutoff)
                likelihood_reads_best_cluster, _, _ = em.expectation(clusters_k[0]['mu'], clusters_k[0]['pi'])
                for c in range(k):
                    reads_best_cluster.loc[:, (k, c + 1)] = likelihood_reads_best_cluster[bitvector.read_inv, c]

            # Save the results
            segs = [path.TopSeg, path.ModSeg, path.SampSeg,
                    path.RefSeg, path.SectSeg, path.ClustMbrSeg]
            fields = dict(top=out_dir, module=path.MOD_CLST,
                          sample=report.sample, ref=report.ref,
                          end5=sect.end5, end3=sect.end3,
                          ext=path.CSVZIP_EXT)
            clust_member_file = path.build(*segs, **fields)
            reads_best_cluster.to_csv(clust_member_file, header=True, index=True,
                                      compression="gzip")
            cluster_out_files.append(clust_member_file)

    return cluster_out_files
