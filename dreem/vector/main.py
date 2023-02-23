import pandas as pd

from dreem.align.reads import BamIndexer
from dreem.util.seq import DNA
from dreem.util.files_sanity import check_library
from dreem.util.cli import *
from dreem.util.path import (BAM_EXT, sanitize, OneRefAlignmentInFilePath,
                             RefsetSeqInFilePath, TopDirPath)
from dreem.vector.mprofile import VectorWriterSpawner


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

    # Change data types.
    bam_paths = list_bam_paths(bamf, bamd)
    coords = list(coords)
    primers = encode_primers(list(primers))

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
                                  parallel=parallel,
                                  max_procs=max_procs,
                                  min_phred=min_phred,
                                  phred_enc=phred_enc,
                                  rerun=rerun)
    return writers.generate_profiles(save_temp, resume)
