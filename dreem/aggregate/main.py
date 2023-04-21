from collections import Counter, defaultdict
from logging import getLogger
from pathlib import Path

import pandas as pd
from click import command

from .library_samples import get_samples_info, get_library_sections
from .mutation_count import process_vectors
from ..util import docdef, path
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
from ..util.sect import encode_primers, get_coords_by_ref
from ..vector.load import VectorLoader

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

    # Extract the arguments
    if library:
        df_library = pd.read_csv(library)
        section_names = get_library_sections(df_library)
        if section_names:
            # Add the coordinates from the library to the coordinates
            # given as a parameter. Remove duplicates and preserve the
            # order by converting the coordinates into dictionary keys
            # first (with Counter) and then casting back to a tuple.
            coords = tuple(Counter(coords + tuple(section_names)))
    else:
        section_names = dict()

    df_samples = check_samples(pd.read_csv(samples)) if samples != "" else None

    # Find all reports among the given report files and directories.
    report_paths = path.find_files_multi(map(Path, mv_file), [path.VecRepSeg])

    # Group coordinates and primers by reference.
    ref_coords = get_coords_by_ref(coords)
    ref_primers = get_coords_by_ref(encode_primers(primers))

    # Aggregate every vector report.
    all_samples = defaultdict(dict)
    for report in report_paths:
        try:
            loader = VectorLoader.load(report)
            sample, ref = loader.sample, loader.ref
            if ref in all_samples[sample]:
                logger.error(f"Got multiple reports for sample '{sample}', "
                             f"ref '{ref}'")
                continue
            report_data = process_vectors(loader,
                                          coords=ref_coords[loader.ref],
                                          primers=ref_primers[loader.ref],
                                          primer_gap=primer_gap,
                                          sect_names=section_names,
                                          out_dir=Path(out_dir))
            all_samples[sample][ref] = report_data
        except Exception as error:
            logger.error(f"Failed to aggregate vectors in {report}: {error}")

    # Add the sample information
    if df_samples is not None:
        for sample in list(all_samples.keys()):
            all_samples[sample].update(get_samples_info(df_samples, sample))
    else:
        for sample in list(all_samples.keys()):
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
