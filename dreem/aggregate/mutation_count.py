from logging import getLogger
from pathlib import Path
from typing import Any

import pandas as pd

from ..calc.auto import quant_muts, KEY_SUB, KEY_HST
from ..vector.load import VectorLoader
from ..util import path
from ..util.sect import Section

logger = getLogger(__name__)


def get_metadata(section: Section):
    return {"section_start": section.end5,
            "section_end": section.end3,
            "sequence": section.seq.decode()}


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
            sect_data[clust] = {KEY_SUB: mus.to_list()}
    return sect_data


def process_vectors(vl: VectorLoader, sections: list[Section], out_dir: Path):
    """ Compute the population average, per-vector, and cluster mutation
    rates (if given) for each section of a set of vectors. Write them to
    CSV files, then return them as a JSON-compatible data structure. """
    # Compute the mutational data for each section.
    per_vect, pop_avgs, clust_mu = quant_muts(vl, sections)
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
            pavg.to_csv(path.build(*segs, **fields, table=path.MUT_POP_AVG))
            # Cluster mutation rates
            if cmus is not None:
                cmus.to_csv(path.build(*segs, **fields, table=path.MUT_CLUSTER))
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
