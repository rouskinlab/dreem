from collections import defaultdict
from logging import getLogger
from pathlib import Path

import pandas as pd
from click import command

from .library_samples import get_samples_info
from .mutation_count import process_vectors
from ..util import docdef
from ..util.cli import (opt_out_dir, opt_temp_dir, opt_save_temp,
                        opt_library, opt_samples,
                        opt_mv_file, opt_clust_file,
                        opt_rnastructure_path, opt_rnastructure_use_temp,
                        opt_rnastructure_fold_args, opt_rnastructure_use_dms,
                        opt_rnastructure_dms_min_unpaired_value,
                        opt_rnastructure_dms_max_paired_value,
                        opt_rnastructure_deltag_ensemble,
                        opt_rnastructure_probability,
                        opt_coords, opt_primers, opt_primer_gap)
from ..util.dependencies import *
from ..util.dump import *
from ..util.files_sanity import check_samples
from ..util.parallel import lock_temp_dir
from ..util.rnastructure import RNAstructure
from ..util.sect import encode_primers
from ..vector.load import open_sections

logger = getLogger(__name__)

params = [
    opt_mv_file,
    opt_coords,
    opt_primers,
    opt_primer_gap,
    opt_library,
    opt_samples,
    opt_clust_file,
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


@command("aggregate", params=params)
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
        min_gap: int,
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
    reports, sections = open_sections(map(Path, mv_file),
                                      coords=coords,
                                      primers=encode_primers(primers),
                                      primer_gap=primer_gap,
                                      library=(Path(library) if library
                                               else None))

    # List the paths for the clustering files.
    clust_paths = list(map(Path, clust_file))

    # Compute the mutation counts for each section of each report.
    all_samples = defaultdict(dict)
    for report in reports:
        try:
            sects_data = process_vectors(report,
                                         sections.list(report.ref),
                                         clust_paths,
                                         Path(out_dir),
                                         min_mut_gap=min_gap)
            all_samples[report.sample][report.ref] = sects_data
        except Exception as error:
            raise
            logger.error(f"Failed to aggregate vectors in {report}: {error}")

    # Add the sample information
    for sample in list(all_samples.keys()):
        if df_samples is not None:
            all_samples[sample].update(get_samples_info(df_samples, sample))
        else:
            all_samples[sample]["sample"] = sample

    rna = RNAstructure(rnastructure_path=rnastructure_path, temp=os.path.join(temp_dir, 'rnastructure'))
    for sample, mut_profiles in all_samples.items():
        for reference in mut_profiles:
            if type(mut_profiles[reference]) is not dict:
                continue

            for section in mut_profiles[reference]:
                if type(mut_profiles[reference][section]) is not dict:
                    continue
                # Add RNAstructure predictions
                mh = rna.run(mut_profiles[reference][section]['sequence'])
                all_samples[sample][reference][section] = {**mut_profiles[reference][section], **mh}
    rna.dump_ledger()

    for sample, mut_profiles in all_samples.items():
        # Write the output
        out = cast_dict(mut_profiles)
        out = sort_dict(out)

        # Make lists in one line
        out = json.dumps(out, cls=NpEncoder, indent=2)
        out = out.replace("]", "[").split("[")
        out = [o.replace("\n  ", "").replace("\n", "") if i % 2 else o for i, o in enumerate(out)]
        out = out[0] + ''.join([('[', ']')[i % 2] + o for i, o in enumerate(out[1:])])

        # Write the output
        with open(os.path.join(out_dir, sample + '.json'), 'w') as f:
            f.write(out)
