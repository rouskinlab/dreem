from functools import cache
from logging import getLogger
from pathlib import Path

import pandas as pd
from click import command

from .library_samples import get_samples_info, get_library_info
from .mutation_count import get_profiles, KEY_SUB
from ..util import docdef, path
from ..util.cli import (opt_out_dir, opt_temp_dir, opt_save_temp,
                        opt_library, opt_samples,
                        opt_bv_files, opt_clustering_file,
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
from ..util.rnastructure import RNAstructure
from ..util.sect import (add_coords_from_library, encode_primers,
                         get_sections, get_coords_by_ref)
from ..util.seq import DNA
from ..vector.load import VectorLoader

logger = getLogger(__name__)

params = [
    opt_bv_files,
    opt_coords,
    opt_primers,
    opt_primer_gap,
    opt_library,
    opt_samples,
    opt_clustering_file,
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


@docdef.auto()
def run(bv_files: tuple[str],
        *,
        out_dir: str,
        temp_dir: str,
        save_temp: bool,
        coords: tuple[tuple[str, int, int], ...],
        primers: tuple[tuple[str, str, str], ...],
        primer_gap: int,
        library: str,
        samples: str,
        clustering_file: str,
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
    # library = check_library(pd.read_csv(library), fasta, out_dir) if library != "" else None
    df_samples = check_samples(pd.read_csv(samples)) if samples != "" else None

    # Find all reports among the given report files and directories.
    report_paths = path.find_files_chain(map(Path, bv_files), [path.VecRepSeg])

    if library:
        # Add coordinates from the library.
        coords = add_coords_from_library(library, coords)

    # Group the coordinates and primers by reference.
    ref_coords = get_coords_by_ref(coords)
    ref_primers = get_coords_by_ref(encode_primers(primers))

    # Ensure consistent reference sequences.
    ref_seqs: dict[str, DNA] = dict()

    # Cache the sections.
    @cache
    def ref_sections(ref: str):
        return get_sections(ref, ref_seqs[ref],
                            coords=ref_coords[ref],
                            primers=ref_primers[ref],
                            primer_gap=primer_gap)

    # Aggregate every vector report.
    all_data = dict()
    for report in report_paths:
        try:
            # Load the report from the file.
            vload = VectorLoader.load(report)
        except Exception as error:
            logger.critical(f"Failed to load vector report {report}: {error}")
            continue
        try:
            if ref_seqs[vload.ref] != vload.seq:
                # The reference sequence does not match another sequence
                # with the same name.
                logger.critical(f"Got >1 sequence for reference '{vload.ref}'")
                continue
        except KeyError:
            # A reference with this name has not yet been encountered.
            ref_seqs[vload.ref] = vload.seq
        # Get the sections of the reference.
        sections = ref_sections(vload.ref)
        # First get information about the reference from the library.
        report_info, sect_names = get_library_info(library, vload.ref)
        # Create a mutational profile for each section.
        try:
            metadata, per_vect, pop_avgs, clusters = get_profiles(vload,
                                                                  sections)  # FIXME add clusters
        except Exception as error:
            logger.critical("Failed to generate mutational profiles for "
                            f"report {report} sections {sections}: {error}")
            continue
        # Add the sample and reference keys to the dict of all samples.
        if vload.sample not in all_data:
            all_data[vload.sample] = dict()
        if vload.ref not in all_data[vload.sample]:
            all_data[vload.sample][vload.ref] = dict()
        # Aggregate the data for every section.
        for sect in sections:
            # Initialize the section's data by merging the metadata and
            # the per-vector data.
            sect_data = metadata[sect.coord] | per_vect[sect.coord]
            # Ensure all metadata are compatible with JSON format.
            for field in list(sect_data.keys()):
                if isinstance(sect_data[field], pd.Series):
                    # Convert Series to dict.
                    sect_data[field] = sect_data[field].to_dict()
            # Add population average data.
            sect_data["pop_avg"] = pop_avgs[sect.coord].to_dict()
            # Add cluster mutation rates.
            for clust, clust_mus in clusters[sect.coord].to_dict().items():
                sect_data[clust] = {KEY_SUB: clust_mus}
            # If the section is named in the library, use that name.
            # Otherwise, use the section's range: end5-end3.
            sect_name = sect_names.get(sect.coord, sect.range)
            # Add the section data to the dict of all information.
            if sect_name in all_data[vload.sample][vload.ref]:
                logger.critical(f"Skipping report file {report} with duplicate "
                                f"sample '{vload.sample}', ref '{vload.ref}', "
                                f"and section '{sect_name}'")
            else:
                all_data[vload.sample][vload.ref][sect_name] = sect_data

    # Add the sample information
    for sample in all_data:
        if df_samples is not None:
            all_data[sample].update(get_samples_info(df_samples, sample))
        else:
            all_data[sample]["sample"] = sample

    rna = RNAstructure(rnastructure_path=rnastructure_path, temp=os.path.join(temp_dir, 'rnastructure'))
    for sample, mut_profiles in all_data.items():
        for reference in mut_profiles:
            if type(mut_profiles[reference]) is not dict:
                continue

            for section in mut_profiles[reference]:
                if type(mut_profiles[reference][section]) is not dict:
                    continue
                # Add RNAstructure predictions
                mh = rna.run(mut_profiles[reference][section]['sequence'])
                all_data[sample][reference][section] = {**mut_profiles[reference][section], **mh}
    rna.dump_ledger()

    for sample, mut_profiles in all_data.items():
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
        print('Outputed to {}'.format(os.path.join(out_dir, sample + '.json')))
    return 1
