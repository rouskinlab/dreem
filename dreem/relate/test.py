from sys import byteorder
import unittest as ut

import numpy as np

from .relate import relate_line
from ..core.rel import as_sam, iter_alignments, NOCOV, MED_QUAL
from ..core.seq import DNA


class TestRelateLineAmbrel(ut.TestCase):
    """ Test function `relate.relate_line`. """

    def assert_equal(self, result: bytearray, expect: np.ndarray):
        with self.subTest(result=result, expect=expect):
            self.assertEqual(result, expect.tobytes())

    @staticmethod
    def relate(ref: str, refseq: DNA, read, qual, cigar, end5):
        """ Generate a SAM line from the given information, and use it
        to compute a relation vector. """
        sam_line = as_sam(b"read", 99, ref, end5, 0, cigar, "=", 1, len(read),
                          read, qual)
        muts = bytearray(NOCOV.to_bytes(1, byteorder) * len(refseq))
        relate_line(sam_line, muts, refseq, len(refseq), ref, MED_QUAL, True)
        return muts

    def iter_cases(self, refseq: DNA):
        """ Iterate through every test case. """
        for read, qual, cigar, end5, end3, relvec in iter_alignments(refseq):
            muts = self.relate("ref", refseq, read, qual, cigar, end5)
            print(refseq, read, qual, cigar, end5, relvec, np.frombuffer(muts, dtype=np.uint8))
            self.assert_equal(muts, relvec)

    def test_a(self):
        self.iter_cases(DNA(b"AA"))


if __name__ == '__main__':
    ut.main()
