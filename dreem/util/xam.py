import itertools
from io import BufferedReader
import os
import re

from dreem.util.cmd import SAMTOOLS_CMD, run_cmd
from dreem.util.dflt import BUFFER_LENGTH
from dreem.util.fa import FastaParser
from dreem.util.ngs import NgsFileBase, SAM_ALIGN_SCORE, SAM_EXTRA_SCORE, SAM_HEADER
from dreem.util.path import OUTPUT_DIR, switch_directory


SAM_EXT = ".sam"
BAM_EXT = ".bam"


class SamBase(NgsFileBase):
    def __init__(self,
                 root_dir: str,
                 ref_file: str,
                 sample: str,
                 xam_file: str) -> None:
        super().__init__(root_dir, ref_file, sample)
        self._ref_names = [ref for ref, _ in FastaParser(ref_file).parse()]
        self._xam_in = xam_file
    
    @property
    def refs(self):
        return self._ref_names
    
    @property
    def xam_in(self):
        return self._xam_in
    
    @property
    def name(self):
        return os.path.splitext(os.path.basename(self.xam_in))[0]
    
    @property
    def sam_name(self):
        return f"{self.name}.sam"

    @property
    def bam_name(self):
        return f"{self.name}{BAM_EXT}"
    
    @property
    def sam_out(self):
        return os.path.join(self.output_dir, self.sam_name)
    
    @property
    def bam_out(self):
        return os.path.join(self.output_dir, self.bam_name)
    
    def _get_bam_split(self, ref: bytes):
        return os.path.join(self.output_dir, f"{ref}{BAM_EXT}")
    
    @staticmethod
    def _index_bam(bam_file: str):
        cmd = [SAMTOOLS_CMD, "index", bam_file]
        run_cmd(cmd)

    def index_bam_in(self):
        self._index_bam(self.xam_in)
    
    def index_bam_out(self):
        self._index_bam(self.bam_out)


class SamRemoveEqualMappers(SamBase):
    _operation_dir = "alignment/3_rem"

    pattern_a = re.compile(SAM_ALIGN_SCORE + rb"(\d+)")
    pattern_x = re.compile(SAM_EXTRA_SCORE + rb"(\d+)")

    def __init__(self, root_dir: str, ref_file: str, sample: str,
                 paired: bool, xam_file: str) -> None:
        super().__init__(root_dir, ref_file, sample, xam_file)
        self._paired = paired
    
    @property
    def paired(self):
        return self._paired

    @staticmethod
    def get_score(line: bytes, ptn: re.Pattern[bytes]):
        return (float(match.groups()[0])
                if (match := ptn.search(line)) else None)
    
    @classmethod
    def is_best_alignment(cls, line: bytes):
        return ((score_x := cls.get_score(line, cls.pattern_x)) is None
                or score_x < cls.get_score(line, cls.pattern_a))
    
    def _iter_paired(self, sam: BufferedReader, line: bytes):
        for line2 in sam:
            if self.is_best_alignment(line) or self.is_best_alignment(line2):
                yield b"".join((line, line2))
            line = sam.readline()
    
    def _iter_single(self, sam: BufferedReader, line: bytes):
        while line:
            if self.is_best_alignment(line):
                yield(line)
            line = sam.readline()

    def _remove_equal_mappers(self, buffer_length=BUFFER_LENGTH):
        with open(self.xam_in, "rb") as sami, open(self.sam_out, "wb") as samo:
            # Copy the header from the input to the output SAM file.
            while (line := sami.readline()).startswith(SAM_HEADER):
                samo.write(line)
            iter_sam = self._iter_paired if self.paired else self._iter_single
            lines = iter_sam(sami, line)
            while text := b"".join(itertools.islice(lines, buffer_length)):
                samo.write(text)
    
    def run(self):
        print("\nRemoving Reads Mapping Equally to Multiple Locations in "
              f"{self.xam_in}\n")
        self._make_output_dir()
        self._remove_equal_mappers()


class SamSorter(SamBase):
    _operation_dir = "alignment/4_sort"

    def _sort(self, name: bool = False):
        cmd = [SAMTOOLS_CMD, "sort"]
        if name:
            cmd.append("-n")
        cmd.extend(["-o", self.bam_out, self.xam_in])
        run_cmd(cmd)
    
    def run(self, name: bool = False):
        print(f"\nSorting {self.xam_in} by Reference and Coordinate\n")
        self._make_output_dir()
        self._sort()


class SamSplitter(SamBase):
    _operation_dir = "alignment/5_split"

    def _get_bam_out_ref(self, ref: bytes):
        return os.path.join(self.output_dir, f"{ref.decode()}{BAM_EXT}")
    
    @property
    def bams_out(self):
        return list(map(self._get_bam_out_ref, self.refs))
    
    @property
    def bam_out(self):
        raise NotImplementedError
    
    def _output_bam_ref(self, ref: bytes, bam_out: str):
        cmd = [SAMTOOLS_CMD, "view", "-b", "-o", bam_out,
               self.xam_in, ref.decode()]
        run_cmd(cmd)
    
    def _split_bam(self):
        self.index_bam_in()
        for ref, bam in zip(self.refs, self.bams_out):
            self._output_bam_ref(ref, bam)
    
    def run(self):
        print(f"\nSplitting {self.xam_in} into Individual References\n")
        self._make_output_dir()
        self._split_bam()


class SamOutputter(SamBase):
    _operation_dir = "alignment"

    @property
    def output_dir(self):
        return self._get_dir(OUTPUT_DIR)
    
    @property
    def bam_out(self):
        return switch_directory(self.xam_in, self.output_dir)
    
    def run(self):
        print(f"\nOutputting Cleaned BAM files to {self.output_dir}\n")
        self._make_output_dir()
        os.rename(self.xam_in, self.bam_out)
