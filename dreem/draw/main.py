

import os

import pandas as pd
from click import command, pass_obj

from ..util import docdef
from ..util.cli import (DreemCommandName, dreem_command,
                        opt_out_dir, opt_draw_input,
                        opt_flat, opt_coords, opt_primers,
                        opt_mutation_fraction, opt_mutation_fraction_identity,
                        opt_base_coverage, opt_mutations_in_barcodes,
                        opt_mutations_per_read_per_sample,  
                        )
from ..util.dump import *


params = [
    opt_draw_input,
    opt_flat,
    opt_coords,
    opt_primers,
    opt_mutation_fraction,
    opt_mutation_fraction_identity,
    opt_base_coverage,
    opt_mutations_in_barcodes,
    opt_mutations_per_read_per_sample,
]


@command(DreemCommandName.DRAW.value, params=params)
# Pass context object.
@pass_obj
# Turn into DREEM command.
@dreem_command()
def cli(**kwargs):
    return run(**kwargs)


@docdef.auto()
def run( 
        input: tuple[str], 
        *,
        flat: list,
        out_dir: str,
        coords: tuple,
        primers: tuple,
        mutation_fraction: bool,
        mutation_fraction_identity:bool,
        base_coverage: bool,
        mutations_in_barcodes: bool,
        mutations_per_read_per_sample: bool,
        ):
    """Run the draw command.

    """




