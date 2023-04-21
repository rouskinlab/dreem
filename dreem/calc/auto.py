from logging import getLogger
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd

from .count import sum_bits, get_bits, QEQ, QEB, QSB
from .mu import cluster_mus
from ..vector.load import VectorLoader
from ..util import path
from ..util.sect import Section
from ..util.seq import (NOCOV, MATCH, DELET, INS_5, INS_3,
                        SUB_A, SUB_C, SUB_G, SUB_T, SUB_N)

logger = getLogger(__name__)

KEY_MAT = "match"
KEY_DEL = "del"
KEY_INS = "ins"
KEY_S_A = "sub_A"
KEY_S_C = "sub_C"
KEY_S_G = "sub_G"
KEY_S_T = "sub_T"
KEY_S_N = "sub_N"
KEY_COV = "cov"
KEY_INF = "info"
KEY_SUB = "sub_rate"
KEY_NAL = "num_aligned"
KEY_MIN = "min_cov"
KEY_HST = "sub_hist"

# Define types of mutations to count by position and by vector.
QUERIES = {
    KEY_MAT: (MATCH | INS_5, QEB),
    KEY_DEL: (DELET, QEQ),
    KEY_INS: (INS_3 | MATCH, QEQ),
    KEY_S_A: (SUB_A, QEQ),
    KEY_S_C: (SUB_C, QEQ),
    KEY_S_G: (SUB_G, QEQ),
    KEY_S_T: (SUB_T, QEQ),
    KEY_S_N: (SUB_N, QEB),
    KEY_COV: (NOCOV, QSB),
}
# Types of mutations to count for each position
Q_BY_POS = list(QUERIES)
# Types of mutations to count for each vector
Q_BY_VEC = [KEY_S_N, KEY_COV]


def get_per_vect(vect_counts: dict[tuple[int, str], pd.Series]):
    # Collect the per-vector information for the section.
    per_vect = dict()
    per_vect[KEY_NAL] = (vect_counts[QUERIES[KEY_COV]] > 0).sum()
    per_vect[KEY_MIN] = vect_counts[QUERIES[KEY_COV]].min()
    hmrg, hmin, hmax = 0.5, 0, vect_counts[QUERIES[KEY_S_N]].max()
    if np.isnan(hmax):
        hmax = 0
    hbins = np.linspace(hmin - hmrg, hmax + hmrg, (hmax - hmin + 1) + 1)
    sub_hist_val = np.histogram(vect_counts[QUERIES[KEY_S_N]], bins=hbins)[0]
    index = np.asarray(np.round((hbins[:-1] + hbins[1:]) / 2), dtype=int)
    per_vect[KEY_HST] = pd.Series(sub_hist_val, index=index)
    return per_vect


def get_pop_avgs(pos_counts: dict[tuple[int, str], pd.Series]):
    pop_avg = {q: pos_counts[QUERIES[q]] for q in Q_BY_POS}
    pop_avg[KEY_INF] = pop_avg[KEY_MAT] + pop_avg[KEY_S_N]
    pop_avg[KEY_SUB] = pop_avg[KEY_S_N] / pop_avg[KEY_INF]
    return pop_avg


def get_clust_mu(loader: VectorLoader, clusters: Path, coord: tuple[int, int]):
    vectors = loader.get_all_vectors(loader.section(*coord).positions)
    return cluster_mus(get_bits(vectors, *QUERIES[KEY_MAT]),
                       get_bits(vectors, *QUERIES[KEY_S_N]),
                       pd.read_csv(clusters))


def get_cluster_coords(loader: VectorLoader, cluster_files: Iterable[Path]):
    """ Find the clustering files with the same sample and reference as
    the vectors and index them by their coordinates. """
    cluster_coords = dict()
    for cfile in cluster_files:
        try:
            fs = path.parse(cfile, path.SampSeg, path.RefSeg,
                            path.SectSeg, path.ClustMbrSeg)
            if fs[path.SAMP] == loader.sample and fs[path.REF] == loader.ref:
                coord = fs[path.END5], fs[path.END3]
                if coord in cluster_coords:
                    logger.warning(f"Skipping duplicate cluster file: {cfile}")
                else:
                    cluster_coords[coord] = cfile
        except Exception as error:
            logger.error(f"Failed to parse cluster file path {cfile}: {error}")
    return cluster_coords


def quant_muts(loader: VectorLoader,
               sections: Iterable[Section],
               cluster_files: Iterable[Path] = ()):
    """
    For every section, gather metadata, count mutations for each vector
    and each position, and compute the cluster mutation rates if one or
    more cluster files are given.

    Parameters
    ----------
    loader: VectorLoader
        Vector loader object
    sections: Iterable[Section]
        One or more sections in which to count the mutations
    cluster_files: Iterable[Path] = ()
        Cluster report files, for computing cluster mus (optional)

    Returns
    -------
    df_row:
        - info: the count per residue of bases for which we have information.
        - cov: the count per residue of bases for which we have coverage.
        - sub_N: the count per residue of bases for which we have mutations.
        - sub_A: the count per residue of bases who are modified to A.
        - sub_C: the count per residue of bases who are modified to C.
        - sub_G: the count per residue of bases who are modified to G.
        - sub_T: the count per residue of bases who are modified to T.
        - sub_rate: the mutation fraction per residue.
        - num_aligned: the number of aligned reads.
    """
    # Count the mutations for each vector and position in each section.
    counts = sum_bits(loader,
                      coords=[sect.coord for sect in sections],
                      by_pos=[QUERIES[q] for q in Q_BY_POS],
                      by_vec=[QUERIES[q] for q in Q_BY_VEC],
                      numeric=True)
    # Get the coordinates of the clusters.
    cluster_coords = get_cluster_coords(loader, cluster_files)
    # Compute mutation rates and other statistics for each section.
    per_vect: dict[tuple[int, int], dict[str, Any]] = dict()
    pop_avgs: dict[tuple[int, int], pd.DataFrame] = dict()
    clust_mu: dict[tuple[int, int], pd.DataFrame] = dict()
    for coord, (vec_counts, pos_counts) in counts.items():
        try:
            # Collect the per-vector information for the section.
            pv = get_per_vect(vec_counts)
            # Collect the population average data for the section.
            pa = pd.DataFrame.from_dict(get_pop_avgs(pos_counts))
            # Compute the mutation rates for each cluster, if any.
            if clust_file := cluster_coords.get(coord):
                clust_mu[coord] = get_clust_mu(loader, clust_file, coord)
        except Exception as error:
            logger.error(f"Failed to compute cluster data for {loader} "
                         f"section {coord}: {error}")
        else:
            # Everything for this section succeeded, so add the data for
            # this section to the collection of data for all sections.
            per_vect[coord] = pv
            pop_avgs[coord] = pa
    if extra_clusts := sorted(set(cluster_coords) - set(clust_mu)):
        logger.warning(f"No sections were given for clusters {extra_clusts}")
    return per_vect, pop_avgs, clust_mu
