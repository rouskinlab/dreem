from logging import getLogger
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd

from ..cluster.files import load_cluster_resps
from ..util import path
from ..util.sect import Section
from ..util.seq import (NOCOV, MATCH, DELET, INS_5, INS_3,
                        SUB_A, SUB_C, SUB_G, SUB_T, SUB_N)
from ..vector.bits import QEQ, QEB, QSB, mvec_to_bvec, sum_bits
from ..vector.load import VectorLoader
from ..vector.mu import cluster_mus

logger = getLogger(__name__)


def get_metadata(section: Section):
    return {"section_start": section.end5,
            "section_end": section.end3,
            "sequence": section.seq.decode()}


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


def summarize_per_vect(vect_counts: dict[tuple[int, str], pd.Series]):
    # Collect the per-vector information for the section.
    per_vect = dict()
    per_vect[KEY_NAL] = (vect_counts[QUERIES[KEY_COV]] > 0).sum()
    per_vect[KEY_MIN] = vect_counts[QUERIES[KEY_COV]].min()
    hmrg, hmin, hmax = 0.5, 0, vect_counts[QUERIES[KEY_S_N]].max()
    if np.isnan(hmax):
        hmax = 0
    hbins = np.linspace(hmin - hmrg, hmax + hmrg, (hmax - hmin + 1) + 1)
    sub_hist_val = np.histogram(vect_counts[QUERIES[KEY_S_N]], bins=hbins)[0]
    index = np.asarray(np.round((hbins[: -1] + hbins[1:]) / 2), dtype=int)
    per_vect[KEY_HST] = pd.Series(sub_hist_val, index=index)
    return per_vect


def summarize_pop_avgs(pos_counts: dict[tuple[int, str], pd.Series]):
    pop_avg = {q: pos_counts[QUERIES[q]] for q in Q_BY_POS}
    pop_avg[KEY_INF] = pop_avg[KEY_MAT] + pop_avg[KEY_S_N]
    pop_avg[KEY_SUB] = pop_avg[KEY_S_N] / pop_avg[KEY_INF]
    return pop_avg


def summarize_clust_mu(loader: VectorLoader,
                       clusters: Path,
                       coord: tuple[int, int],
                       min_gap: int):
    vectors = loader.get_all_vectors(loader.section(*coord).positions)
    return cluster_mus(mvec_to_bvec(vectors, *QUERIES[KEY_MAT]),
                       mvec_to_bvec(vectors, *QUERIES[KEY_S_N]),
                       load_cluster_resps(clusters),
                       min_gap)


def get_cluster_coords(loader: VectorLoader, cluster_files: Iterable[Path]):
    """ Find the clustering files with the same sample and reference as
    the vectors and index them by their coordinates. """
    cluster_coords = dict()
    for cfile in cluster_files:
        try:
            fs = path.parse(cfile, path.SampSeg, path.RefSeg,
                            path.SectSeg, path.ClustTabSeg)
            if fs[path.SAMP] == loader.sample and fs[path.REF] == loader.ref:
                coord = fs[path.END5], fs[path.END3]
                if coord in cluster_coords:
                    logger.warning(f"Skipping duplicate cluster file: {cfile}")
                else:
                    cluster_coords[coord] = cfile
        except Exception as error:
            logger.error(f"Failed to parse cluster file path {cfile}: {error}")
    return cluster_coords


def summarize(loader: VectorLoader,
              sections: Iterable[Section],
              cluster_files: Iterable[Path] = (),
              min_gap: int | None = None):
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
    min_gap: int | None = None
        Minimum separation between two mutations; must be an integer ≥ 0
        if cluster_files is not empty, otherwise ignored

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
    if cluster_coords and min_gap is None:
        raise TypeError("min_gap must be given if any cluster_files are given")
    # Compute mutation rates and other statistics for each section.
    per_vect: dict[tuple[int, int], dict[str, Any]] = dict()
    pop_avgs: dict[tuple[int, int], pd.DataFrame] = dict()
    clust_mu: dict[tuple[int, int], pd.DataFrame] = dict()
    for coord, (vec_counts, pos_counts) in counts.items():
        # Collect the per-vector information for the section.
        pvec = summarize_per_vect(vec_counts)
        # Collect the population average data for the section.
        pavg = pd.DataFrame.from_dict(summarize_pop_avgs(pos_counts))
        try:
            # Compute the mutation rates for each cluster, if any.
            if cfile := cluster_coords.get(coord):
                clust_mu[coord] = summarize_clust_mu(loader, cfile, coord, min_gap)
        except Exception as error:
            raise
            logger.error(f"Failed to compute cluster data for {loader} "
                         f"section {coord}: {error}")
        else:
            # Everything for this section succeeded, so add the data for
            # this section to the collection of data for all sections.
            per_vect[coord] = pvec
            pop_avgs[coord] = pavg
    if extra_clusts := sorted(set(cluster_coords) - set(clust_mu)):
        logger.warning(f"No sections were given for clusters: {extra_clusts}")
    return per_vect, pop_avgs, clust_mu


def jsonify_section(metadata: dict[str, Any],
                    per_vect: dict[str, Any],
                    pop_avgs: pd.DataFrame,
                    clust_mu: pd.DataFrame | None):
    """ Convert the metadata and mutation data for a section from arrays
    to dictionaries that can be saved in JSON format. """
    # Initialize the section's data by merging the metadata and
    # the per-vector data.
    sect_data = metadata | per_vect
    # Ensure all metadata are compatible with JSON format.
    for field in list(sect_data.keys()):
        if isinstance(sect_data[field], pd.Series):
            # Convert Series to dict.
            sect_data[field] = sect_data[field].to_list()
    # Add population average data.
    sect_data["pop_avg"] = {field: values.to_list()
                            for field, values in pop_avgs.items()}
    # Add cluster mutation rates.
    if clust_mu is not None:
        for clust, mus in clust_mu.to_dict().items():
            sect_data[clust] = {KEY_SUB: mus}
    return sect_data


def process_vectors(vl: VectorLoader,
                    sections: list[Section],
                    cluster_files: list[Path],
                    out_dir: Path,
                    min_mut_gap: int | None = None):
    """ Compute the population average, per-vector, and cluster mutation
    rates (if given) for each section of a set of vectors. Write them to
    CSV files, then return them as a JSON-compatible data structure. """
    # Compute the mutational data for each section.
    per_vect, pop_avgs, clust_mu = summarize(vl, sections, cluster_files,
                                             min_mut_gap)
    # JSON-ify the data for every section.
    json_data = dict()
    for sect in sections:
        # Get the mutation data for the section.
        meta = get_metadata(sect)
        pvec = per_vect[sect.coord]
        pavg = pop_avgs[sect.coord]
        cmus = clust_mu.get(sect.coord)
        # Write the mutation data to CSV files.
        try:
            segs = [path.ModSeg, path.SampSeg, path.RefSeg, path.SectSeg]
            fields = {path.TOP: out_dir,
                      path.MOD: path.MOD_AGGR,
                      path.SAMP: vl.sample,
                      path.REF: vl.ref,
                      path.END5: sect.end5,
                      path.END3: sect.end3}
            # Make the parent directory, if it does not exist.
            path.build(*segs, **fields).mkdir(parents=True, exist_ok=True)
            # Create the CSV files.
            segs.append(path.MutTabSeg)
            fields[path.EXT] = path.CSV_EXT
            # Histogram of mutations per vector
            pvec[KEY_HST].to_csv(path.build(*segs, **fields,
                                            table=path.MUT_PER_VEC))
            # Population average reactivities
            pavg.to_csv(path.build(*segs, **fields,
                                   table=path.MUT_POP_AVG))
            # Cluster mutation rates
            if cmus is not None:
                cmus.to_csv(path.build(*segs, **fields,
                                       table=path.MUT_CLUSTER))
        except Exception as error:
            logger.error(f"Failed to write mutation data for {sect}: {error}")
        # Convert the data to a JSON-compatible data structure.
        if sect.name in json_data:
            logger.warning(f"Skipping duplicate section: {sect}")
            continue
        try:
            json_data[sect.name] = jsonify_section(meta, pvec, pavg, cmus)
        except Exception as error:
            logger.error(f"Failed to make {sect} JSON-compatible: {error}")
    return json_data
