from logging import getLogger
from pathlib import Path

import pandas as pd
from click import command

from ..call.load import CallVecLoader
from ..cluster.load import ClusterLoader
from ..core import docdef, path
from ..core.cli import (opt_out_dir, opt_temp_dir, opt_save_temp,
                        opt_library, opt_samples,
                        opt_mvec, opt_clust,
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
    opt_mvec,
    opt_coords,
    opt_primers,
    opt_primer_gap,
    opt_library,
    opt_samples,
    opt_clust,
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


@command(path.MOD_TABLE, params=params)
def cli(**kwargs):
    return run(**kwargs)


@docdef.auto()
def run(call_report: tuple[str, ...],
        clust_report: tuple[str, ...]):
    """
    Run the table module.
    """
    for report in clust_report
