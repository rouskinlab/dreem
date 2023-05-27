from logging import getLogger
from pathlib import Path

import pandas as pd
from click import command

from ..core import docdef, path
from ..core.cli import (opt_temp_dir, opt_save_temp, opt_table,
                        opt_rnastructure_use_temp,
                        opt_rnastructure_fold_args, opt_rnastructure_use_dms,
                        opt_rnastructure_dms_min_unpaired_value,
                        opt_rnastructure_dms_max_paired_value,
                        opt_rnastructure_deltag_ensemble,
                        opt_rnastructure_probability)
from ..core.dependencies import *
from ..core.dump import *
from ..core.files_sanity import check_samples
from ..core.parallel import lock_temp_dir
from ..core.rna import RNAstructure
from ..core.sect import encode_primers
from ..relate.load import open_sections as open_vectors

logger = getLogger(__name__)

params = [
    opt_table,
    opt_rnastructure_use_temp,
    opt_rnastructure_fold_args,
    opt_rnastructure_use_dms,
    opt_rnastructure_dms_min_unpaired_value,
    opt_rnastructure_dms_max_paired_value,
    opt_rnastructure_deltag_ensemble,
    opt_rnastructure_probability,
    opt_temp_dir,
    opt_save_temp,
]


@command(path.MOD_FOLD, params=params)
def cli(**kwargs):
    return run(**kwargs)


@lock_temp_dir
@docdef.auto()
def run(table: tuple[str, ...],
        *,
        temp_dir: str,
        save_temp: bool,
        rnastructure_use_temp: bool,
        rnastructure_fold_args: str,
        rnastructure_use_dms: str,
        rnastructure_dms_min_unpaired_value: float,
        rnastructure_dms_max_paired_value: float,
        rnastructure_deltag_ensemble: bool,
        rnastructure_probability: bool):
    """
    Run the fold module.
    """

    check_rnastructure_exists()
