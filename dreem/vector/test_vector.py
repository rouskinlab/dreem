import os
import itertools
import unittest
from unittest import TestCase

from dreem.util.util import *
from vector import *



class TestConsensusByte(TestCase):
    pass



class TestSamFlag(TestCase):
    flags = ("paired", "proper", "unmap", "munmap", "rev", "mrev",
             "first", "second", "secondary", "qcfail", "dup", "supp")

    bits = len(flags)

    # attributes

    def test_flag_attrs(self):
        sf = SamFlag(0)
        self.assertListEqual(sf.__slots__, list(self.flags))
        self.assertEqual(sf.MAX_FLAG, 4095)

    # valid flags

    def test_flag_0(self):
        sf = SamFlag(0)
        self.assertFalse(any(sf.__getattribute__(flag) for flag in self.flags))

    def test_flag_99(self):
        sf = SamFlag(99)
        self.assertTrue(all([sf.paired, sf.proper, sf.mrev, sf.first]))
        self.assertFalse(any([sf.unmap, sf.munmap, sf.rev, sf.second,
                              sf.secondary, sf.qcfail, sf.dup, sf.supp]))

    def test_flag_4095(self):
        sf = SamFlag(4095)
        self.assertTrue(all(sf.__getattribute__(flag) for flag in self.flags))

    def test_flags_all_valid(self):
        for flag, flags in enumerate(itertools.product([False, True],
                                                       repeat=self.bits)):
            sf = SamFlag(flag)
            sf_flags = list(map(sf.__getattribute__, self.flags))
            real_flags = list(reversed(flags))
            self.assertListEqual(sf_flags, real_flags)

    # invalid flags

    def test_flag_4096(self):
        with self.assertRaises(ValueError):
            SamFlag(4096)

    def test_flag_neg1(self):
        with self.assertRaises(ValueError):
            SamFlag(-1)

    def test_flag_float(self):
        with self.assertRaises(TypeError):
            SamFlag(99.0)

    def test_flag_str(self):
        with self.assertRaises(TypeError):
            SamFlag("99")

    def test_flag_bytes(self):
        with self.assertRaises(TypeError):
            SamFlag(b"99")


class TestSamRead(TestCase):
    # valid SAM lines

    def test_read_end5(self):
        line = b"FS10000136:97:BPN80019-0831:1:1101:6870:1040	99	end5	1	22	31M65S	=	1	-96	CAGCACTCAGAGCTAATACGACTCACTATAGATAATTGTGTACAAAGTAGAGATGTATCCAATTATGTGACTACCTTTGTGTAATAAAAATTTGTT	,FFFFFF:FF,:FFFFF:::,FF,F,FFF,:F:F,FF:F,F,:F,::F,FF,:,:FFF:F,FFF:FFF:,F,F::FF::FFF,F:FF::,FFF,:F	MD:Z:31	XG:i:0	NM:i:0	XM:i:0	XN:i:0	XO:i:0	AS:i:62	YS:i:62	YT:Z:CP"
        sr = SamRead(line)
        sf99 = SamFlag(99)
        self.assertEqual(sr.qname, b"FS10000136:97:BPN80019-0831:1:1101:6870:1040")
        self.assertListEqual(
            list(map(sr.flag.__getattribute__, TestSamFlag.flags)),
            list(map(sf99.__getattribute__, TestSamFlag.flags)))
        self.assertEqual(sr.rname, b"end5")
        self.assertEqual(sr.pos, 1)
        self.assertEqual(sr.cigar, b"31M65S")
        self.assertEqual(sr.tlen, -96)
        self.assertEqual(sr.seq, b"CAGCACTCAGAGCTAATACGACTCACTATAGATAATTGTGTACAAAGTAGAGATGTATCCAATTATGTGACTACCTTTGTGTAATAAAAATTTGTT")
        self.assertEqual(len(sr), 96)

    # invalid SAM lines

    def test_seq_qual_different_lengths(self):
        with self.assertRaises(ValueError):
            SamRead(b"Q	0	R	1	100	4=	*	*	4	ACGT	III")

    def test_read_empty(self):
        with self.assertRaises(ValueError):
            SamRead(b"")

    def test_read_str(self):
        with self.assertRaises(TypeError):
            line = "FS10000136:97:BPN80019-0831:1:1101:6870:1040	99	end5	1	22	31M65S	=	1	-96	CAGCACTCAGAGCTAATACGACTCACTATAGATAATTGTGTACAAAGTAGAGATGTATCCAATTATGTGACTACCTTTGTGTAATAAAAATTTGTT	,FFFFFF:FF,:FFFFF:::,FF,F,FFF,:F:F,FF:F,F,:F,::F,FF,:,:FFF:F,FFF:FFF:,F,F::FF::FFF,F:FF::,FFF,:F	MD:Z:31	XG:i:0	NM:i:0	XM:i:0	XN:i:0	XO:i:0	AS:i:62	YS:i:62	YT:Z:CP"
            SamRead(line)


class TestParseCigar(TestCase):
    # valid CIGAR strings

    def test_cigar_1op(self):
        parsed = list(parse_cigar(b"300="))
        self.assertListEqual(parsed, [(b"=", 300)])

    def test_cigar_5ops(self):
        parsed = list(parse_cigar(b"5I1X65=1D8S"))
        self.assertListEqual(parsed,
            [(b"I", 5), (b"X", 1), (b"=", 65), (b"D", 1), (b"S", 8)])

    # invalid CIGAR strings

    def test_cigar_wrong_order(self):
        with self.assertRaises(ValueError):
            list(parse_cigar(b"=300"))

    def test_cigar_double_op(self):
        with self.assertRaises(ValueError):
            list(parse_cigar(b"300=="))

    def test_cigar_invalid_op(self):
        with self.assertRaises(ValueError):
            list(parse_cigar(b"300x"))

    def test_cigar_missing_op(self):
        with self.assertRaises(ValueError):
            list(parse_cigar(b"300"))

    def test_cigar_neg_length(self):
        with self.assertRaises(ValueError):
            list(parse_cigar(b"-300="))

    def test_cigar_missing_length(self):
        with self.assertRaises(ValueError):
            list(parse_cigar(b"="))

    def test_cigar_zero_length(self):
        with self.assertRaises(ValueError):
            list(parse_cigar(b"0="))


class TestCompMutsRead(TestCase):
    def test_valid_span_match(self):
        ref = b"ACGT"
        first, last = 1, 4
        line = b"Q	0	R	1	100	4=	*	*	4	ACGT	IIII"
        expect = MATCH*4
        muts = vectorize_read(ref, first, last, SamRead(line))
        self.assertTrue(muts == expect)
    
    def test_valid_span_1del(self):
        ref = b"ACGT"
        first, last = 1, 4
        line = b"Q	0	R	1	100	1=1D2=	*	*	4	AGT	III"
        expect = MATCH + DELET + MATCH*2
        muts = vectorize_read(ref, first, last, SamRead(line))
        self.assertTrue(muts == expect)

    def test_valid_span_1del_ambig(self):
        ref = b"ACCCCCT"
        first, last = 1, 7
        line = b"Q	0	R	1	100	3=1D3=	*	*	4	ACCCCT	IIIIII"
        expect = MATCH + MADEL*5 + MATCH
        muts = vectorize_read(ref, first, last, SamRead(line))
        self.assertTrue(muts == expect)

    def test_valid_span_2del_ambig(self):
        ref = b"ACCGGT"
        first, last = 1, 6
        line = b"Q	0	R	1	100	2=2D2=	*	*	4	ACGT	IIII"
        expect = MATCH + MADEL*4 + MATCH
        muts = vectorize_read(ref, first, last, SamRead(line))
        self.assertTrue(muts == expect)

    def test_valid_span_2del_ambig_2(self):
        ref = b"ACCGGT"
        first, last = 1, 6
        line = b"Q	0	R	1	100	2=1D3=	*	*	4	ACGGT	IIIII"
        expect = MATCH + MADEL*2 + MATCH*3
        muts = vectorize_read(ref, first, last, SamRead(line))
        self.assertTrue(muts == expect)
    
    def test_valid_span_tunnel2(self):
        ref = b"GATATG"
        first, last = 1, 6
        line = b"Q	0	R	1	100	2=2D2=	*	*	4	GATG	IIII"
        expect = MATCH + MADEL*4 + MATCH
        muts = vectorize_read(ref, first, last, SamRead(line))
        self.assertTrue(muts == expect)


class TestCompMutsPairedReads(TestCase):
    def test_valid_no_muts(self):
        ref = b"ACGT"
        first, last = 1, 4
        line1 = b"Q	0	R	1	30	4=	*	*	4	ACGT	IIII"
        line2 = line1
        expect = MATCH * 4
        muts = vectorize_pair(ref, first, last, SamRead(line1), SamRead(line2))
        self.assertTrue(muts == expect)


if __name__ == "__main__":
    unittest.main()
