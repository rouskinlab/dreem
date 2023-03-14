import sys
import os

from ..demultiplex.demultiplex import demultiplex_run


def run(library_csv,demulti_workspace,mixed_fastq1,mixed_fastq2,clipped:int=0,index_tolerance:int=0,parallel:bool=False,mismatch_tolerence:int=0):
    x=demultiplex_run( library_csv=library_csv,
                                 demulti_workspace=demulti_workspace,
                                 mixed_fastq1=mixed_fastq1,
                                 mixed_fastq2=mixed_fastq2,
                                 clipped=clipped,
                                 index_tolerance=index_tolerance,
                                 parallel=parallel,
                                 mismatch_tolerence=mismatch_tolerence)
    return x


if __name__ == '__main__':
    pass
