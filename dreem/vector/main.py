import pandas as pd

from ..align.reads import BamIndexer
from ..util.cli import *
from ..util.files_sanity import check_library
from ..util.path import (BAM_EXT, sanitize, OneRefAlignmentInFilePath,
                         RefsetSeqInFilePath, TopDirPath)
from ..util.seq import DNA
from ..vector.mprofile import VectorWriterSpawner


def add_coords_from_library(library_path: str,
                            coords: list[tuple[str, int, int]],
                            fasta: str,
                            out_dir: str):
    library = check_library(pd.read_csv(library_path), fasta, out_dir)
    for construct, end5, end3 in zip(library["construct"],
                                     library["section_start"],
                                     library["section_end"], strict=True):
        coord = (str(os.path.splitext(construct)[0]), int(end5), int(end3))
        if coord not in coords:
            coords.append(coord)


def encode_primers(primers: list[tuple[str, str, str]]):
    return [(ref, DNA(fwd.encode()), DNA(rev.encode()))
            for ref, fwd, rev in primers]


def list_bam_paths(bamf: tuple[str, ...], bamd: tuple[str, ...]):
    bam_paths: dict[str, OneRefAlignmentInFilePath] = dict()

    def add_bam_path(bam):
        bam_path = sanitize(bam)
        if not os.path.isfile(bam_path):
            logging.error(f"Skipping non-existant BAM file: {bam_path}")
            return
        if not bam_path.endswith(BAM_EXT):
            logging.error(f"Skipping non-BAM-formatted file: {bam_path}")
            return
        if bam_path in bam_paths:
            logging.warning(f"Skipping duplicate BAM file: {bam_path}")
            return
        bam_paths[bam_path] = OneRefAlignmentInFilePath.parse(bam_path)

    for bam_file in bamf:
        add_bam_path(bam_file)
    for bam_dir in bamd:
        if not os.path.isdir(bam_dir):
            logging.error(f"No such BAM directory: {bam_dir}")
            continue
        for bam_file in os.listdir(bam_dir):
            if bam_file.endswith(BAM_EXT):
                add_bam_path(os.path.join(bam_dir, bam_file))

    return list(bam_paths.values())


@click.command(DreemCommandName.VECTOR.value)
# Input files
@opt_fasta
@opt_bamf
@opt_bamd
# SAM options
@opt_phred_enc
@opt_min_phred
# Output directories
@opt_out_dir
@opt_temp_dir
# File generation
@opt_rerun
@opt_resume
@opt_save_temp
# Parallelization
@opt_parallel
@opt_max_procs
@opt_batch_size
# Regions
@opt_library
@opt_cfill
@opt_coords
@opt_primers
@opt_primer_gap
# Pass context object
@click.pass_obj
# Turn into DREEM command
@dreem_command(imports=("fasta", "bamf"),
               exports="mp_report")
def run(out_dir: str,
        temp_dir: str,
        bamf: tuple[str],
        bamd: tuple[str],
        fasta: str,
        coords: tuple[tuple[str, int, int], ...],
        primers: tuple[tuple[str, str, str], ...],
        primer_gap: int,
        cfill: bool,
        library: str,
        batch_size: float,
        parallel: bool,
        max_procs: int,
        phred_enc: int,
        min_phred: int,
        rerun: bool,
        save_temp: bool,
        resume: bool):
    """
    Run the vectoring step.
    Generate a vector encoding mutations for each read
    (or read pair, if paired-end).

    FASTA (path): reference sequence(s) in FASTA format

    BAM_FILES (paths): list of one or more alignment files (space separated)
                       in BAM format
    """

    # List out all the paths to BAM files.
    bam_paths = list_bam_paths(bamf, bamd)
    # Convert coords to list (expected type for VectorWriterSpawner).
    coords = list(coords)
    # Convert given primer sequences (str) to DNA objects.
    primers = encode_primers(list(primers))
    # Convert batch_size from mebibytes (2^20 bytes) to bytes.
    bytes_per_batch = round(batch_size * 2 ** 20)

    # Index every BAM file.
    for bam_path in bam_paths:
        BamIndexer(xam=bam_path, num_cpus=max_procs, resume=resume).run()

    # If a library file is given, add coordinates from the file.
    if library:
        add_coords_from_library(library_path=library,
                                coords=coords,
                                fasta=fasta,
                                out_dir=out_dir)

    # Compute mutation vectors for each BAM file.
    writers = VectorWriterSpawner(out_dir=TopDirPath.parse(out_dir),
                                  temp_dir=TopDirPath.parse(temp_dir),
                                  refset_path=RefsetSeqInFilePath.parse(fasta),
                                  bam_paths=bam_paths,
                                  coords=coords,
                                  primers=primers,
                                  primer_gap=primer_gap,
                                  cfill=cfill,
                                  batch_size=bytes_per_batch,
                                  parallel=parallel,
                                  max_procs=max_procs,
                                  min_phred=min_phred,
                                  phred_enc=phred_enc,
                                  rerun=rerun)
    return writers.generate_profiles(save_temp, resume)
