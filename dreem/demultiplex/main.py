import sys
import os

from ..demultiplex.demultiplex import demultiplex_run

from click import command, pass_obj

from ..util.cli import (
    opt_barcode_length, opt_barcode_start, opt_parallel_demultiplexing, opt_clipped_demultiplexing,
    opt_mismatch_tolerence, opt_index_tolerence, opt_demulti_overwrite, opt_fasta, opt_library, opt_fastq1, opt_fastq2)

params = [
    # Inputs
    opt_fasta,
    opt_fastq1,
    opt_fastq2,
    opt_library,
    opt_barcode_start,
    opt_barcode_length,
    

    # options
    opt_parallel_demultiplexing,
    opt_clipped_demultiplexing,
    opt_mismatch_tolerence,
    opt_index_tolerence,
    opt_demulti_overwrite

]


# Turn into DREEM command.

@command("demultiplex", params=params)
# Pass context object.
@pass_obj
def cli(**kwargs):
    return run(**kwargs)


def run(library_csv: str, demulti_workspace: str, mixed_fastq1: str, mixed_fastq2: str, fasta: str, barcode_start=0,
        barcode_length=0, clipped: int = 0, index_tolerance: int = 0, parallel: bool = False,
        mismatch_tolerence: int = 0, overwrite: bool = False, secondary_signiture_filtering:bool=False):
    return demultiplex_run(library_csv=library_csv,
                           overwrite=overwrite,
                           demulti_workspace=demulti_workspace,
                           mixed_fastq1=mixed_fastq1,
                           mixed_fastq2=mixed_fastq2,
                           barcode_start=barcode_start,
                           barcode_length=barcode_length,
                           clipped=clipped,
                           index_tolerance=index_tolerance,
                           parallel=parallel,
                           fasta=fasta,
                           mismatch_tolerence=mismatch_tolerence,
                           secondary_signiture_filtering=secondary_signiture_filtering)


if __name__ == '__main__':
    pass
