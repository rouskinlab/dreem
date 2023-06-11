from itertools import product
from sys import byteorder
import unittest as ut

from .relate import relate_line
from ..core.rel import as_sam, iter_alignments, NOCOV, MED_QUAL
from ..core.seq import DNA, BASES


class TestRelateLineAmbrel(ut.TestCase):
    """ Test function `relate.relate_line`. """

    @staticmethod
    def relate(ref: str, refseq: DNA, read, qual, cigar, end5):
        """ Generate a SAM line from the given information, and use it
        to compute a relation vector. """
        sam_line = as_sam(b"read", 99, ref, end5, 0, cigar, "=", 1, len(read),
                          read, qual)
        muts = bytearray(NOCOV.to_bytes(1, byteorder) * len(refseq))
        relate_line(sam_line, muts, refseq, len(refseq), ref, MED_QUAL, True)
        return muts

    def iter_cases(self, refseq: DNA, max_ins: int = 2):
        """ Iterate through every test case. """
        for read, qual, cigar, end5, end3, relvec in iter_alignments(refseq,
                                                                     max_ins,
                                                                     max_ins,
                                                                     max_ins):
            result = self.relate("ref", refseq, read, qual, cigar, end5)
            with self.subTest(relvec=relvec, result=result):
                self.assertEqual(relvec, result)

    def test_aaaa_1ins(self):
        """ Test all possible reads with ≤ 1 insertion from AAAA. """
        self.iter_cases(DNA(b"AAAA"), 1)

    def test_aacc_1ins(self):
        """ Test all possible reads with ≤ 1 insertion from AACC. """
        self.iter_cases(DNA(b"AACC"), 1)

    def test_acgt_1ins(self):
        """ Test all possible reads with ≤ 1 insertion from ACGT. """
        self.iter_cases(DNA(b"ACGT"), 1)



if __name__ == '__main__':
    ut.main()
