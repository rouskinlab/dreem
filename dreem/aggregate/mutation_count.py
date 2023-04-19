from logging import getLogger
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd

from ..vector.calc import cluster_mus, sum_bits, get_bits, QEQ, QEB, QEP, QSB
from ..vector.load import VectorLoader
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
    KEY_INS: (INS_3, QEP),
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


def get_metadata(section: Section):
    return {"section_start": section.end5,
            "section_end": section.end3,
            "sequence": section.seq.decode()}


def get_per_vect(vect_counts: dict[tuple[int, str], pd.Series]):
    # Collect the per-vector information for the section.
    per_vect = dict()
    per_vect[KEY_NAL] = (vect_counts[QUERIES[KEY_COV]] > 0).sum()
    per_vect[KEY_MIN] = vect_counts[QUERIES[KEY_COV]].min()
    hmrg, hmin, hmax = 0.5, 0, vect_counts[QUERIES[KEY_S_N]].max()
    if np.isnan(hmax):
        hmax = 0
    hbins = np.linspace(hmin - hmrg, hmax + hmrg, (hmax - hmin + 1) + 1)
    sub_hist_val, sub_hist_idx = np.histogram(vect_counts[QUERIES[KEY_S_N]],
                                              bins=hbins)
    per_vect[KEY_HST] = pd.Series(sub_hist_val, index=sub_hist_idx[:-1])
    return per_vect


def get_pop_avgs(pos_counts: dict[tuple[int, str], pd.Series]):
    pop_avg = {q: pos_counts[QUERIES[q]] for q in Q_BY_POS}
    pop_avg[KEY_INF] = pop_avg[KEY_MAT] + pop_avg[KEY_S_N]
    pop_avg[KEY_SUB] = pop_avg[KEY_S_N] / pop_avg[KEY_INF]
    return pop_avg


def get_clust_mu(loader: VectorLoader, clust_file: Path):
    vectors = loader.get_all_vectors(section.positions)
    return cluster_mus(get_bits(vectors, *QUERIES[KEY_MAT]),
                       get_bits(vectors, *QUERIES[KEY_S_N]),
                       pd.read_csv(clust_file))


def get_profiles(loader: VectorLoader,
                 coords: Iterable[tuple[int, int]],
                 cluster_files: Iterable[Path] = ()):
    """
    Generate a mutational profile from mutation vectors and sections.

    Parameters
    ----------
    loader: VectorLoader
        Vector loader object
    coords: Iterable[tuple[int, int]]
        One or more pairs of 5' and 3' coordinates between which to
        count the number of mutations.
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
    if cluster_files is None:
        cluster_files = list()
    cluster_coords = {}  # FIXME: add parsing of cluster report paths
    # Count the mutations.
    counts = sum_bits(loader, coords=coords,
                      by_pos=[QUERIES[q] for q in Q_BY_POS],
                      by_vec=[QUERIES[q] for q in Q_BY_VEC],
                      numeric=True)
    # Construct the mutational profiles.
    metadata: dict[tuple[int, int], dict[str, Any]] = dict()
    per_vect: dict[tuple[int, int], dict[str, Any]] = dict()
    pop_avgs: dict[tuple[int, int], pd.DataFrame] = dict()
    clust_mu: dict[tuple[int, int], pd.DataFrame] = dict()
    for coord, (vec_counts, pos_counts) in counts.items():
        # Collect the general information for the section.
        metadata[coord] = get_metadata(loader.section(*coord))
        # Collect the per-vector information for the section.
        per_vect[coord] = get_per_vect(vec_counts)
        # Collect the population average data for the section.
        pop_avgs[coord] = pd.DataFrame.from_dict(get_pop_avgs(pos_counts))
        # Compute the mutation rates for each cluster of the section.
        if clust_file := cluster_coords.get(coord):
            try:
                clust_mu[coord] = get_clust_mu(loader, clust_file)
            except Exception as error:
                logger.error(f"Failed to compute cluster data for {loader} "
                             f"section {coord}: {error}")
    # FIXME if extra_clusts := { cluster_coords}
    return metadata, per_vect, pop_avgs, clust_mu
