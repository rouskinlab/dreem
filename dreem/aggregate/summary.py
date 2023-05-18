from logging import getLogger
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from ..cluster.load import ClusterLoader
from ..util import path
from ..util.sect import Section
from ..util.seq import (NOCOV, MATCH, DELET, INS_5, INS_3,
                        SUB_A, SUB_C, SUB_G, SUB_T, SUB_N)
from ..vector.bits import sum_bits, QEQ, QEB, QSB
from ..vector.load import VectorLoader

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


def jsonify_section(metadata: dict[str, Any],
                    per_vect: dict[str, Any],
                    pop_avgs: pd.DataFrame):
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
    return sect_data


def vectors(loader: VectorLoader, sections: list[Section], out_dir: Path):
    """ Compute the population average and per-vector mutation rates for
    each section of a set of vectors. Write them to CSV files, and then
    return them as a JSON-compatible data structure. """
    # Count the mutations for each vector and position in each section.
    counts = sum_bits(loader,
                      coords=[sect.coord for sect in sections],
                      by_pos=[QUERIES[q] for q in Q_BY_POS],
                      by_vec=[QUERIES[q] for q in Q_BY_VEC],
                      numeric=True)
    # Compute mutation rates and other statistics for each section.
    per_vect: dict[tuple[int, int], dict[str, Any]] = dict()
    pop_avgs: dict[tuple[int, int], pd.DataFrame] = dict()
    for coord, (vec_counts, pos_counts) in counts.items():
        # Collect the per-vector information for the section.
        per_vect[coord] = summarize_per_vect(vec_counts)
        # Collect the population average data for the section.
        pop_avgs[coord] = pd.DataFrame.from_dict(summarize_pop_avgs(pos_counts))
    # JSON-ify the data for every section.
    json_data = dict()
    for section in sections:
        # Get the mutation data for the section.
        meta = get_metadata(section)
        pvec = per_vect[section.coord]
        pavg = pop_avgs[section.coord]
        # Write the mutation data to CSV files.
        try:
            segs = [path.ModSeg, path.SampSeg, path.RefSeg, path.SectSeg]
            fields = {path.TOP: out_dir,
                      path.MOD: path.MOD_AGGR,
                      path.SAMP: loader.sample,
                      path.REF: loader.ref,
                      path.SECT: section.name}
            # Make the parent directory, if it does not exist.
            path.build(*segs, **fields).mkdir(parents=True, exist_ok=True)
            # Create the CSV files.
            segs.append(path.MutTabSeg)
            fields[path.EXT] = path.CSV_EXT
            # Histogram of mutations per vector.
            pvec[KEY_HST].to_csv(path.build(*segs, **fields,
                                            table=path.MUT_PER_VEC))
            # Population average reactivities.
            pavg.to_csv(path.build(*segs, **fields,
                                   table=path.MUT_POP_AVG))
        except Exception as error:
            logger.error(f"Failed to write mutation data for {section}: {error}")
        # Convert the data to a JSON-compatible data structure.
        if section.name in json_data:
            logger.warning(f"Skipping duplicate section: {section}")
            continue
        try:
            json_data[section.name] = jsonify_section(meta, pvec, pavg)
        except Exception as error:
            logger.error(f"Failed to make {section} JSON-compatible: {error}")
    return json_data


def clusters(loader: ClusterLoader, out_dir: Path):
    # Write the mutation rates and proportions to CSV files.
    try:
        segs = [path.ModSeg, path.SampSeg, path.RefSeg, path.SectSeg,
                path.ClustTabSeg]
        fields = {path.TOP: out_dir,
                  path.MOD: path.MOD_AGGR,
                  path.SAMP: loader.sample,
                  path.REF: loader.ref,
                  path.SECT: loader.sect,
                  path.RUN: 0,
                  path.EXT: path.CSV_EXT}
        # Mutation rates.
        loader.mus.to_csv(path.buildpar(*segs, **fields,
                                        table=path.CLUST_MUS_TABLE))
        # Cluster proportions.
        loader.props.to_csv(path.buildpar(*segs, **fields,
                                          table=path.CLUST_PROP_TABLE))
    except Exception as error:
        logger.error(f"Failed to write cluster mus for {loader}: {error}")
    # Get the metadata for the section.
    json_data = get_metadata(loader.section)
    # Add the cluster mutation rates and proportions.
    if not loader.mus.columns.equals(loader.props.index):
        raise ValueError(f"Indexes of mutation rates and proportions disagree")
    for k, i in loader.mus.columns:
        json_data[f"K{k}-{i}"] = {"mus": loader.mus[k, i].to_dict(),
                                  "prop": loader.props[k, i]}
    return json_data
