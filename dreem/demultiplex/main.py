from click import command

from ..demultiplex.demultiplex import demultiplex_run
from ..util.cli import (
    opt_barcode_length, opt_barcode_start, opt_parallel_demultiplexing, opt_clipped_demultiplexing,
    opt_mismatch_tolerence, opt_index_tolerence, opt_demulti_overwrite, opt_fasta, opt_library, opt_fastq1, opt_fastq2,opt_out_dir,opt_temp_dir)

params = [
    # Inputs
    opt_fasta,
    opt_fastq1,
    opt_fastq2,
    opt_library,
    opt_barcode_start,
    opt_barcode_length,
    opt_out_dir,
    opt_temp_dir,

    # options
    opt_parallel_demultiplexing,
    opt_clipped_demultiplexing,
    opt_mismatch_tolerence,
    opt_index_tolerence,
    opt_demulti_overwrite,


]


# Turn into DREEM command.

@command("demultiplex", params=params)
def cli(**kwargs):
    return run(**kwargs)


def run(library: str, out_dir: str, temp_dir:str, fastq1: str, fastq2: str, fasta: str, barcode_start=0,
        barcode_length=0, clipped: int = 0, index_tolerance: int = 0, parallel_demultiplexing: bool = False,
        mismatch_tolerence: int = 0, demulti_overwrite: bool = False):
    
    return demultiplex_run(library_csv=library,
                           overwrite=demulti_overwrite,
                           demulti_workspace=temp_dir,
                           report_folder=out_dir,
                           mixed_fastq1=fastq1,
                           mixed_fastq2=fastq2,
                           barcode_start=barcode_start,
                           barcode_length=barcode_length,
                           clipped=clipped,
                           index_tolerance=index_tolerance,
                           parallel=parallel_demultiplexing,
                           fasta=fasta,
                           mismatch_tolerence=mismatch_tolerence)


if __name__ == '__main__':
    pass
