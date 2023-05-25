from logging import getLogger
from pathlib import Path

import pandas as pd
from click import command

from .library_samples import get_samples_info
from .summary import clusters, vectors
from ..cluster.load import ClusterLoader
from ..core import docdef, path
from ..core.cli import (opt_out_dir, opt_temp_dir, opt_save_temp,
                        opt_samples, opt_table,
                        opt_rnastructure_path, opt_rnastructure_use_temp,
                        opt_rnastructure_fold_args, opt_rnastructure_use_dms,
                        opt_rnastructure_dms_min_unpaired_value,
                        opt_rnastructure_dms_max_paired_value,
                        opt_rnastructure_deltag_ensemble,
                        opt_rnastructure_probability,
                        opt_coords, opt_primers, opt_primer_gap)
from ..core.dependencies import *
from ..core.dump import *
from ..core.files_sanity import check_samples
from ..core.parallel import lock_temp_dir
from ..core.rnastructure import RNAstructure
from ..core.sect import encode_primers
from ..relate.load import open_sections as open_vectors

logger = getLogger(__name__)

params = [
    opt_table,
    opt_samples,
    opt_rnastructure_path,
    opt_rnastructure_use_temp,
    opt_rnastructure_fold_args,
    opt_rnastructure_use_dms,
    opt_rnastructure_dms_min_unpaired_value,
    opt_rnastructure_dms_max_paired_value,
    opt_rnastructure_deltag_ensemble,
    opt_rnastructure_probability,
    opt_out_dir,
    opt_temp_dir,
    opt_save_temp,
]


@command(path.MOD_AGGR, params=params)
def cli(**kwargs):
    return run(**kwargs)


@lock_temp_dir
@docdef.auto()
def run(mv_file: tuple[str, ...],
        clust_file: tuple[str, ...],
        *,
        out_dir: str,
        temp_dir: str,
        save_temp: bool,
        coords: tuple[tuple[str, int, int], ...],
        primers: tuple[tuple[str, str, str], ...],
        primer_gap: int,
        library: str,
        samples: str,
        rnastructure_path: str,
        rnastructure_use_temp: bool,
        rnastructure_fold_args: str,
        rnastructure_use_dms: str,
        rnastructure_dms_min_unpaired_value: float,
        rnastructure_dms_max_paired_value: float,
        rnastructure_deltag_ensemble: bool,
        rnastructure_probability: bool):
    """
    Run the aggregate module.
    """

    check_rnastructure_exists(rnastructure_path)

    df_samples = check_samples(pd.read_csv(samples)) if samples != "" else None

    # Open all vector reports and get the sections for each.
    loaders, sects = open_vectors(map(Path, mv_file),
                                  coords=coords,
                                  primers=encode_primers(primers),
                                  primer_gap=primer_gap,
                                  library=(Path(library) if library
                                           else None))
    # Summarize each section of each vectoring report.
    summary = dict()
    for ld in loaders:
        try:
            if ld.sample not in summary:
                summary[ld.sample] = dict()
            summary[ld.sample][ld.ref] = vectors(ld,
                                                 sects.list(ld.ref),
                                                 Path(out_dir))
        except Exception as error:
            raise
            logger.error(f"Failed to aggregate vectors in {ld}: {error}")

    # Summarize each clustering run.
    for cluster_report in clust_file:
        try:
            ld = ClusterLoader.open(Path(cluster_report))
            if ld.sample not in summary:
                summary[ld.sample] = dict()
            if ld.ref not in summary[ld.sample]:
                summary[ld.sample][ld.ref] = dict()
            if ld.sect not in summary[ld.sample][ld.ref]:
                summary[ld.sample][ld.ref][ld.sect] = dict()
            summary[ld.sample][ld.ref][ld.sect].update(clusters(ld, Path(out_dir)))
        except Exception as error:
            raise
            logger.error(
                f"Failed to aggregate clusters in {cluster_report}: {error}")

    # Add the sample information
    for sample in list(summary.keys()):
        if df_samples is not None:
            summary[sample].update(get_samples_info(df_samples, sample))
        else:
            summary[sample]["sample"] = sample

    rna = RNAstructure(rnastructure_path=rnastructure_path, temp=os.path.join(temp_dir, 'rnastructure'))
    for sample, mut_profiles in summary.items():
        for reference in mut_profiles:
            if type(mut_profiles[reference]) is not dict:
                continue

            for section in mut_profiles[reference]:
                if type(mut_profiles[reference][section]) is not dict:
                    continue
                # Add RNAstructure predictions
                mh = rna.run(mut_profiles[reference][section]['sequence'])
                summary[sample][reference][section] = {**mut_profiles[reference][section], **mh}
    rna.dump_ledger()

    for sample, mut_profiles in summary.items():
        # Write the output
        out = cast_dict(mut_profiles)
        #out = sort_dict(out)

        # Make lists in one line
        out = json.dumps(out, cls=NpEncoder, indent=2)
        out = out.replace("]", "[").split("[")
        out = [o.replace("\n  ", "").replace("\n", "") if i % 2 else o for i, o in enumerate(out)]
        out = out[0] + ''.join([('[', ']')[i % 2] + o for i, o in enumerate(out[1:])])

        # Write the output
        with open(os.path.join(out_dir, sample + '.json'), 'w') as f:
            f.write(out)
