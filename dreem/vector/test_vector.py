import os
import itertools
import unittest
from unittest import TestCase

from dreem.util.util import *
from dreem.vector.vector import *



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
        read = SamRead(line)
        sf99 = SamFlag(99)
        self.assertEqual(read.qname, b"FS10000136:97:BPN80019-0831:1:1101:6870:1040")
        self.assertListEqual(
            list(map(read.flag.__getattribute__, TestSamFlag.flags)),
            list(map(sf99.__getattribute__, TestSamFlag.flags)))
        self.assertEqual(read.rname, b"end5")
        self.assertEqual(read.pos, 1)
        self.assertEqual(read.cigar, b"31M65S")
        self.assertEqual(read.seq, b"CAGCACTCAGAGCTAATACGACTCACTATAGATAATTGTGTACAAAGTAGAGATGTATCCAATTATGTGACTACCTTTGTGTAATAAAAATTTGTT")
        self.assertEqual(len(read), 96)

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


class TestVectorizeReadMatchSub(TestCase):
    def test_valid_span_match_eq(self):
        ref = b"ACGT"
        first, last = 1, 4
        line = b"Q	0	R	1	100	4=	*	*	4	ACGT	IIII"
        expect = MATCH*4
        muts = vectorize_read(ref, first, last, SamRead(line))
        self.assertTrue(muts == expect)
    
    def test_valid_span_match_m(self):
        ref = b"ACGT"
        first, last = 1, 4
        line = b"Q	0	R	1	100	4M	*	*	4	ACGT	IIII"
        expect = MATCH*4
        muts = vectorize_read(ref, first, last, SamRead(line))
        self.assertTrue(muts == expect)


class TestVectorizeReadOneDel(TestCase):
    def test_1del_5prm_match_match(self):
        ref = b"ACCG"
        first, last = 1, 4
        line = b"Q	0	R	1	100	1M1D2M	*	*	4	ACG	III"
        expect = MATCH + MADEL*2 + MATCH
        muts = vectorize_read(ref, first, last, SamRead(line))
        self.assertTrue(muts == expect)

    def test_1del_3prm_match_match(self):
        ref = b"TGGC"
        first, last = 1, 4
        line = b"Q	0	R	1	100	2M1D1M	*	*	4	TGC	III"
        expect = MATCH + MADEL*2 + MATCH
        muts = vectorize_read(ref, first, last, SamRead(line))
        self.assertTrue(muts == expect)

    def test_1del_5prm_match_match_low_qual(self):
        ref = b"ACCG"
        first, last = 1, 4
        line = b"Q	0	R	1	100	1M1D2M	*	*	4	ACG	I!I"
        DSM_C = (MADEL[0]|(SUB_N[0]^SUB_C[0])).to_bytes()
        expect = MATCH + DSM_C*2 + MATCH
        muts = vectorize_read(ref, first, last, SamRead(line))
        self.assertTrue(muts == expect)

    def test_1del_3prm_match_match_low_qual(self):
        ref = b"TGGC"
        first, last = 1, 4
        line = b"Q	0	R	1	100	2M1D1M	*	*	4	TGC	I!I"
        DSM_G = (MADEL[0]|(SUB_N[0]^SUB_G[0])).to_bytes()
        expect = MATCH + DSM_G*2 + MATCH
        muts = vectorize_read(ref, first, last, SamRead(line))
        self.assertTrue(muts == expect)

    def test_1del_5prm_mismatch_mismatch_same(self):
        ref = b"ACCG"
        first, last = 1, 4
        line = b"Q	0	R	1	100	1M1D2M	*	*	4	ATG	III"
        DEL_T = (DELET[0]|SUB_T[0]).to_bytes()
        expect = MATCH + DEL_T*2 + MATCH
        muts = vectorize_read(ref, first, last, SamRead(line))
        self.assertTrue(muts == expect)

    def test_1del_3prm_mismatch_mismatch_same(self):
        ref = b"TGGC"
        first, last = 1, 4
        line = b"Q	0	R	1	100	2M1D1M	*	*	4	TAC	III"
        DEL_A = (DELET[0]|SUB_A[0]).to_bytes()
        expect = MATCH + DEL_A*2 + MATCH
        muts = vectorize_read(ref, first, last, SamRead(line))
        self.assertTrue(muts == expect)

    def test_1del_5prm_mismatch_mismatch_diff(self):
        ref = b"ACGT"
        first, last = 1, 4
        line = b"Q	0	R	1	100	1M1D2M	*	*	4	ATT	III"
        DEL_T = (DELET[0]|SUB_T[0]).to_bytes()
        expect = MATCH + DEL_T*2 + MATCH
        muts = vectorize_read(ref, first, last, SamRead(line))
        self.assertTrue(muts == expect)

    def test_1del_3prm_mismatch_mismatch_diff(self):
        ref = b"TGCA"
        first, last = 1, 4
        line = b"Q	0	R	1	100	2M1D1M	*	*	4	TAA	III"
        DEL_A = (DELET[0]|SUB_A[0]).to_bytes()
        expect = MATCH + DEL_A*2 + MATCH
        muts = vectorize_read(ref, first, last, SamRead(line))
        self.assertTrue(muts == expect)

    def test_1del_5prm_match_mismatch(self):
        ref = b"ACGT"
        first, last = 1, 4
        line = b"Q	0	R	1	100	1M1D2M	*	*	4	AGT	III"
        expect = MATCH + DELET + MATCH*2
        muts = vectorize_read(ref, first, last, SamRead(line))
        self.assertTrue(muts == expect)

    def test_1del_3prm_match_mismatch(self):
        ref = b"TCGA"
        first, last = 1, 4
        line = b"Q	0	R	1	100	2M1D1M	*	*	4	TCA	III"
        expect = MATCH*2 + DELET + MATCH
        muts = vectorize_read(ref, first, last, SamRead(line))
        self.assertTrue(muts == expect)

    def test_1del_5prm_match_mismatch_low_qual(self):
        ref = b"ACGT"
        first, last = 1, 4
        line = b"Q	0	R	1	100	1M1D2M	*	*	4	AGT	I!I"
        DSM_C = (MADEL[0]|(SUB_N[0]^SUB_C[0])).to_bytes()
        DSM_G = (MADEL[0]|(SUB_N[0]^SUB_G[0])).to_bytes()
        expect = MATCH + DSM_C + DSM_G + MATCH
        muts = vectorize_read(ref, first, last, SamRead(line))
        self.assertTrue(muts == expect)

    def test_1del_3prm_match_mismatch_low_qual(self):
        ref = b"CATG"
        first, last = 1, 4
        line = b"Q	0	R	1	100	2M1D1M	*	*	4	CAG	I!I"
        DSM_A = (MADEL[0]|(SUB_N[0]^SUB_A[0])).to_bytes()
        DSM_T = (MADEL[0]|(SUB_N[0]^SUB_T[0])).to_bytes()
        expect = MATCH + DSM_A + DSM_T + MATCH
        muts = vectorize_read(ref, first, last, SamRead(line))
        self.assertTrue(muts == expect)

    def test_1del_5prm_mismatch_match(self):
        ref = b"ACGT"
        first, last = 1, 4
        line = b"Q	0	R	1	100	1M1D2M	*	*	4	ACT	III"
        expect = MATCH + DELET + SUB_C + MATCH
        muts = vectorize_read(ref, first, last, SamRead(line))
        self.assertTrue(muts == expect)

    def test_1del_3prm_mismatch_match(self):
        ref = b"ACGT"
        first, last = 1, 4
        line = b"Q	0	R	1	100	2M1D1M	*	*	4	AGT	III"
        expect = MATCH + SUB_G + DELET + MATCH
        muts = vectorize_read(ref, first, last, SamRead(line))
        self.assertTrue(muts == expect)


class TestVectorizeReadMultiDels(TestCase):
    def test_2of2dels(self):
        ref = b"ACCG"
        first, last = 1, 4
        line = b"Q	0	R	1	100	1M2D1M	*	*	4	AG	II"
        expect = MATCH + DELET*2 + MATCH
        muts = vectorize_read(ref, first, last, SamRead(line))
        self.assertTrue(muts == expect)

    def test_2of3dels_1(self):
        ref = b"ACCCG"
        first, last = 1, 5
        line = b"Q	0	R	1	100	1M2D2M	*	*	4	ACG	III"
        expect = MATCH + MADEL*3 + MATCH
        muts = vectorize_read(ref, first, last, SamRead(line))
        self.assertTrue(muts == expect)

    def test_2of3dels_2(self):
        ref = b"ACCCG"
        first, last = 1, 5
        line = b"Q	0	R	1	100	1M1D1M1D1M	*	*	4	ACG	III"
        expect = MATCH + MADEL*3 + MATCH
        muts = vectorize_read(ref, first, last, SamRead(line))
        self.assertTrue(muts == expect)

    def test_2of3dels_3(self):
        ref = b"ACCCG"
        first, last = 1, 5
        line = b"Q	0	R	1	100	2M2D1M	*	*	4	ACG	III"
        expect = MATCH + MADEL*3 + MATCH
        muts = vectorize_read(ref, first, last, SamRead(line))
        self.assertTrue(muts == expect)

    def test_2of4dels_1(self):
        ref = b"ACCGGT"
        first, last = 1, 6
        line = b"Q	0	R	1	100	1M1D2M1D1M	*	*	4	ACGT	IIII"
        expect = MATCH + MADEL*4 + MATCH
        muts = vectorize_read(ref, first, last, SamRead(line))
        self.assertTrue(muts == expect)

    def test_2of4dels_2(self):
        ref = b"ACCGGT"
        first, last = 1, 6
        line = b"Q	0	R	1	100	2M2D2M	*	*	4	ACGT	IIII"
        expect = MATCH + MADEL*4 + MATCH
        muts = vectorize_read(ref, first, last, SamRead(line))
        self.assertTrue(muts == expect)

    def test_2of6dels(self):
        ref = b"ACCCCCCG"
        first, last = 1, 8
        line = b"Q	0	R	1	100	3M2D3M	*	*	4	ACCCCG	IIIIII"
        expect = MATCH + MADEL*6 + MATCH
        muts = vectorize_read(ref, first, last, SamRead(line))
        self.assertTrue(muts == expect)
    
    def test_tunnel_2_5prm(self):
        ref = b"ACGCGCGT"
        first, last = 1, 8
        line = b"Q	0	R	1	100	1M2D5M	*	*	4	ACGCGT	IIIIII"
        expect = MATCH + MADEL*6 + MATCH
        muts = vectorize_read(ref, first, last, SamRead(line))
        self.assertTrue(muts == expect)

    def test_tunnel_2_mid(self):
        ref = b"ACGCGCGT"
        first, last = 1, 8
        line = b"Q	0	R	1	100	3M2D3M	*	*	4	ACGCGT	IIIIII"
        expect = MATCH + MADEL*6 + MATCH
        muts = vectorize_read(ref, first, last, SamRead(line))
        self.assertTrue(muts == expect)

    def test_tunnel_2_3prm(self):
        ref = b"ACGCGCGT"
        first, last = 1, 8
        line = b"Q	0	R	1	100	5M2D1M	*	*	4	ACGCGT	IIIIII"
        expect = MATCH + MADEL*6 + MATCH
        muts = vectorize_read(ref, first, last, SamRead(line))
        self.assertTrue(muts == expect)

    def test_tunnel_3_5prm(self):
        ref = b"TCAGCAGCAT"
        first, last = 1, 10
        line = b"Q	0	R	1	100	1M3D6M	*	*	4	TCAGCAT	IIIIIII"
        expect = MATCH + MADEL*8 + MATCH
        muts = vectorize_read(ref, first, last, SamRead(line))
        self.assertTrue(muts == expect)

    def test_tunnel_3_mid(self):
        ref = b"TCAGCAGCAT"
        first, last = 1, 10
        line = b"Q	0	R	1	100	2M3D5M	*	*	4	TCAGCAT	IIIIIII"
        expect = MATCH + MADEL*8 + MATCH
        muts = vectorize_read(ref, first, last, SamRead(line))
        self.assertTrue(muts == expect)

    def test_tunnel_3_3prm(self):
        ref = b"TCAGCAGCAT"
        first, last = 1, 10
        line = b"Q	0	R	1	100	6M3D1M	*	*	4	TCAGCAT	IIIIIII"
        expect = MATCH + MADEL*8 + MATCH
        muts = vectorize_read(ref, first, last, SamRead(line))
        self.assertTrue(muts == expect)


class TestVectorizeReadOneIns(TestCase):
    def test_1ins(self):
        ref = b"ACTT"
        first, last = 1, 4
        line = b"Q	0	R	1	100	2M1I2M	*	*	4	ACGTT	IIIII"
        MIN_5 = (MATCH[0]|INS_5[0]).to_bytes()
        MIN_3 = (MATCH[0]|INS_3[0]).to_bytes()
        expect = MATCH + MIN_5 + MIN_3 + MATCH
        muts = vectorize_read(ref, first, last, SamRead(line))
        self.assertTrue(muts == expect)

    def test_1ins_5prm_match_match(self):
        ref = b"ACT"
        first, last = 1, 3
        line = b"Q	0	R	1	100	1M1I2M	*	*	4	ACCT	IIII"
        MIN_5 = (MATCH[0]|INS_5[0]).to_bytes()
        MIN_3 = (MATCH[0]|INS_3[0]).to_bytes()
        MIN_B = (MIN_5[0]|MIN_3[0]).to_bytes()
        expect = MIN_5 + MIN_B + MIN_3
        muts = vectorize_read(ref, first, last, SamRead(line))
        self.assertTrue(muts == expect)

    def test_1ins_3prm_match_match(self):
        ref = b"ACT"
        first, last = 1, 3
        line = b"Q	0	R	1	100	2M1I1M	*	*	4	ACCT	IIII"
        MIN_5 = (MATCH[0]|INS_5[0]).to_bytes()
        MIN_3 = (MATCH[0]|INS_3[0]).to_bytes()
        MIN_B = (MIN_5[0]|MIN_3[0]).to_bytes()
        expect = MIN_5 + MIN_B + MIN_3
        muts = vectorize_read(ref, first, last, SamRead(line))
        self.assertTrue(muts == expect)
    
    def test_1ins_5prm_match_match_match(self):
        ref = b"ACCT"
        first, last = 1, 4
        line = b"Q	0	R	1	100	1M1I3M	*	*	5	ACCCT	IIIII"
        MIN_5 = (MATCH[0]|INS_5[0]).to_bytes()
        MIN_3 = (MATCH[0]|INS_3[0]).to_bytes()
        MIN_B = (MIN_5[0]|MIN_3[0]).to_bytes()
        expect = MIN_5 + MIN_B*2 + MIN_3
        muts = vectorize_read(ref, first, last, SamRead(line))
        self.assertTrue(muts == expect)

    def test_1ins_center_match_match_match(self):
        ref = b"ACCT"
        first, last = 1, 4
        line = b"Q	0	R	1	100	2M1I2M	*	*	5	ACCCT	IIIII"
        MIN_5 = (MATCH[0]|INS_5[0]).to_bytes()
        MIN_3 = (MATCH[0]|INS_3[0]).to_bytes()
        MIN_B = (MIN_5[0]|MIN_3[0]).to_bytes()
        expect = MIN_5 + MIN_B*2 + MIN_3
        muts = vectorize_read(ref, first, last, SamRead(line))
        self.assertTrue(muts == expect)

    def test_1ins_3prm_match_match_match(self):
        ref = b"ACCT"
        first, last = 1, 4
        line = b"Q	0	R	1	100	3M1I1M	*	*	5	ACCCT	IIIII"
        MIN_5 = (MATCH[0]|INS_5[0]).to_bytes()
        MIN_3 = (MATCH[0]|INS_3[0]).to_bytes()
        MIN_B = (MIN_5[0]|MIN_3[0]).to_bytes()
        expect = MIN_5 + MIN_B*2 + MIN_3
        muts = vectorize_read(ref, first, last, SamRead(line))
        self.assertTrue(muts == expect)

    def test_1ins_5prm_match_lowq(self):
        ref = b"ACTA"
        first, last = 1, 4
        line = b"Q	0	R	1	100	1M1I3M	*	*	4	ACGTA	II!II"
        POS_1 = (MATCH[0]|INS_5[0]).to_bytes()
        POS_2 = (MATCH[0]|INS_5[0]|INS_3[0]|(ANY_N[0]^SUB_C[0])).to_bytes()
        POS_3 = (MATCH[0]|INS_3[0]|INS_5[0]|(ANY_N[0]^SUB_T[0])).to_bytes()
        POS_4 = (MATCH[0]|INS_3[0]).to_bytes()
        expect = POS_1 + POS_2 + POS_3 + POS_4
        muts = vectorize_read(ref, first, last, SamRead(line))
        self.assertTrue(muts == expect)

    def test_1ins_mid_match_lowq(self):
        ref = b"ACTA"
        first, last = 1, 4
        line = b"Q	0	R	1	100	2M1I2M	*	*	4	ACGTA	II!II"
        POS_1 = (MATCH[0]|INS_5[0]).to_bytes()
        POS_2 = (MATCH[0]|INS_5[0]|INS_3[0]|(ANY_N[0]^SUB_C[0])).to_bytes()
        POS_3 = (MATCH[0]|INS_3[0]|INS_5[0]|(ANY_N[0]^SUB_T[0])).to_bytes()
        POS_4 = (MATCH[0]|INS_3[0]).to_bytes()
        expect = POS_1 + POS_2 + POS_3 + POS_4
        muts = vectorize_read(ref, first, last, SamRead(line))
        self.assertTrue(muts == expect)

    def test_1ins_3prm_match_lowq(self):
        ref = b"ACTA"
        first, last = 1, 4
        line = b"Q	0	R	1	100	3M1I1M	*	*	4	ACGTA	II!II"
        POS_1 = (MATCH[0]|INS_5[0]).to_bytes()
        POS_2 = (MATCH[0]|INS_5[0]|INS_3[0]|(ANY_N[0]^SUB_C[0])).to_bytes()
        POS_3 = (MATCH[0]|INS_3[0]|INS_5[0]|(ANY_N[0]^SUB_T[0])).to_bytes()
        POS_4 = (MATCH[0]|INS_3[0]).to_bytes()
        expect = POS_1 + POS_2 + POS_3 + POS_4
        muts = vectorize_read(ref, first, last, SamRead(line))
        self.assertTrue(muts == expect)

    def test_1ins_5prm_lowq_lowq(self):
        """
        Possible alignments:
        Ref:    A-CATA    AC-ATA    ACA-TA    ACAT-A
        Read:   ACNNTA    ACNNTA    ACNNTA    ACNNTA
        """
        ref = b"ACATA"
        first, last = 1, 5
        line = b"Q	0	R	1	100	2M1I3M	*	*	6	ACGATA	II!!II"
        POS_1 = (MATCH[0]|INS_5[0]).to_bytes()
        POS_2 = (MATCH[0]|INS_5[0]|INS_3[0]|(ANY_N[0]^SUB_C[0])).to_bytes()
        POS_3 = (MATCH[0]|INS_5[0]|INS_3[0]|(ANY_N[0]^SUB_A[0])).to_bytes()
        POS_4 = (MATCH[0]|INS_5[0]|INS_3[0]|(ANY_N[0]^SUB_T[0])).to_bytes()
        POS_5 = (MATCH[0]|INS_3[0]).to_bytes()
        expect = POS_1 + POS_2 + POS_3 + POS_4 + POS_5
        muts = vectorize_read(ref, first, last, SamRead(line))
        self.assertTrue(muts == expect)
    
    def test_1ins_3prm_lowq_lowq(self):
        """
        Possible alignments:
        Ref:    A-CGTA    AC-GTA    ACG-TA    ACGT-A
        Read:   ACNNTA    ACNNTA    ACNNTA    ACNNTA
        """
        ref = b"ACGTA"
        first, last = 1, 5
        line = b"Q	0	R	1	100	2M1I3M	*	*	6	ACGATA	II!!II"
        POS_1 = (MATCH[0]|INS_5[0]).to_bytes()
        POS_2 = (MATCH[0]|INS_5[0]|INS_3[0]|(ANY_N[0]^SUB_C[0])).to_bytes()
        POS_3 = (MATCH[0]|INS_5[0]|INS_3[0]|(ANY_N[0]^SUB_G[0])).to_bytes()
        POS_4 = (MATCH[0]|INS_5[0]|INS_3[0]|(ANY_N[0]^SUB_T[0])).to_bytes()
        POS_5 = (MATCH[0]|INS_3[0]).to_bytes()
        expect = POS_1 + POS_2 + POS_3 + POS_4 + POS_5
        muts = vectorize_read(ref, first, last, SamRead(line))
        self.assertTrue(muts == expect)

    def test_1ins_5prm_mismatch_match(self):
        """
        Possible alignments:
        Read:   CAGT
        Ref:    C-AT
        """
        ref = b"CAT"
        first, last = 1, 3
        line = b"Q	0	R	1	100	1M1I2M	*	*	4	CAGT	IIII"
        POS_1 = (MATCH[0]|INS_5[0]).to_bytes()
        POS_2 = (SUB_G[0]|INS_3[0]).to_bytes()
        POS_3 = MATCH
        expect = POS_1 + POS_2 + POS_3
        muts = vectorize_read(ref, first, last, SamRead(line))
        self.assertTrue(muts == expect)
    
    def test_1ins_3prm_mismatch_match(self):
        """
        Possible alignments:
        Read:   CAGT
        Ref:    CG-T
        """
        ref = b"CGT"
        first, last = 1, 3
        line = b"Q	0	R	1	100	2M1I1M	*	*	4	CAGT	IIII"
        POS_1 = MATCH
        POS_2 = (SUB_A[0]|INS_5[0]).to_bytes()
        POS_3 = (MATCH[0]|INS_3[0]).to_bytes()
        expect = POS_1 + POS_2 + POS_3
        muts = vectorize_read(ref, first, last, SamRead(line))
        self.assertTrue(muts == expect)

    def test_1ins_5prm_mismatch_mismatch(self):
        """
        Possible alignments:
        Read:   CATC    CATC
        Ref:    C-GC    CG-C
        """
        ref = b"CGC"
        first, last = 1, 3
        line = b"Q	0	R	1	100	1M1I2M	*	*	4	CATC	IIII"
        POS_1 = (MATCH[0]|INS_5[0]).to_bytes()
        POS_2 = (SUB_A[0]|SUB_T[0]|INS_5[0]|INS_3[0]).to_bytes()
        POS_3 = (MATCH[0]|INS_3[0]).to_bytes()
        expect = POS_1 + POS_2 + POS_3
        muts = vectorize_read(ref, first, last, SamRead(line))
        self.assertTrue(muts == expect)

    def test_1ins_3prm_mismatch_mismatch(self):
        """
        Possible alignments:
        Read:   CATC    CATC
        Ref:    C-GC    CG-C
        """
        ref = b"CGC"
        first, last = 1, 3
        line = b"Q	0	R	1	100	2M1I1M	*	*	4	CATC	IIII"
        POS_1 = (MATCH[0]|INS_5[0]).to_bytes()
        POS_2 = (SUB_A[0]|SUB_T[0]|INS_5[0]|INS_3[0]).to_bytes()
        POS_3 = (MATCH[0]|INS_3[0]).to_bytes()
        expect = POS_1 + POS_2 + POS_3
        muts = vectorize_read(ref, first, last, SamRead(line))
        self.assertTrue(muts == expect)

    def test_1ins_centered(self):
        """
        Possible alignments:
        Read:   AAATAAA
        Ref:    AAA-AAA
        """
        ref = b"AAAAAA"
        first, last = 1, 6
        line = b"Q	0	R	1	100	3M1I3M	*	*	4	AAATAAA	IIIIIII"
        POS_3 = (MATCH[0]|INS_5[0]).to_bytes()
        POS_4 = (MATCH[0]|INS_3[0]).to_bytes()
        expect = MATCH*2 + POS_3 + POS_4 + MATCH*2
        muts = vectorize_read(ref, first, last, SamRead(line))
        self.assertTrue(muts == expect)

    def test_1ins_off_centered_5prm(self):
        """
        Possible alignments:
        Read:   AAATAAA    AAATAAA    AAATAAA
        Ref:    -AAAAAA    A-AAAAA    AA-AAAA
        """
        ref = b"AAAAAA"
        first, last = 1, 6
        line = b"Q	0	R	1	100	1M1I5M	*	*	4	AAATAAA	IIIIIII"
        MI_53 = (MATCH[0]|INS_5[0]|INS_3[0]).to_bytes()
        POS_3 = (SUB_T[0]|INS_3[0]).to_bytes()
        expect = MI_53*2 + POS_3 + MATCH*3
        muts = vectorize_read(ref, first, last, SamRead(line))
        self.assertTrue(muts == expect)

    def test_1ins_off_centered_3prm(self):
        """
        Possible alignments:
        Read:   AAATAAA    AAATAAA    AAATAAA
        Ref:    AAAA-AA    AAAAA-A    AAAAAA-
        """
        ref = b"AAAAAA"
        first, last = 1, 6
        line = b"Q	0	R	1	100	5M1I1M	*	*	4	AAATAAA	IIIIIII"
        MI_53 = (MATCH[0]|INS_5[0]|INS_3[0]).to_bytes()
        POS_4 = (SUB_T[0]|INS_5[0]).to_bytes()
        expect = MATCH*3 + POS_4 + MI_53*2
        muts = vectorize_read(ref, first, last, SamRead(line))
        self.assertTrue(muts == expect)


class TestVectorizeReadMultiInns(TestCase):
    def test_2ins_indep(self):
        """
        Possible alignments:
        Read:   ATCGAT
        Ref:    A-CG-T
        """
        ref = b"ACGT"
        first, last = 1, 4
        line = b"Q	0	R	1	100	1M1I2M1I1M	*	*	6	ATCGAT	IIIIII"
        MIN_5 = (MATCH[0]|INS_5[0]).to_bytes()
        MIN_3 = (MATCH[0]|INS_3[0]).to_bytes()
        expect = MIN_5 + MIN_3 + MIN_5 + MIN_3
        muts = vectorize_read(ref, first, last, SamRead(line))
        self.assertTrue(muts == expect)

    def test_2ins_5prm_ymdt230103(self):
        ref = b"GAGGATCTGCCTGTATCGTCCGGAG"
        first, last = 1, 25
        line = b"Q	0	R	1	100	7M1I14M1I4M	*	*	27	GAGGATCTTGCCTGTATCGTCCGGGAG	IIIIIIIIIIIIIIIIIIIIIIIIIII"
        MIN_5 = (MATCH[0]|INS_5[0]).to_bytes()
        MIN_3 = (MATCH[0]|INS_3[0]).to_bytes()
        MIN_B = (MIN_5[0]|MIN_3[0]).to_bytes()
        expect = MATCH*6 + MIN_5 + MIN_B + MIN_3 + MATCH*11 + MIN_5 + MIN_B*2 + MIN_3 + MATCH
        muts = vectorize_read(ref, first, last, SamRead(line))
        self.assertTrue(muts == expect)

    def test_2ins_3prm_ymdt230103(self):
        ref = b"GAGGATCTGCCTGTATCGTCCGGAG"
        first, last = 1, 25
        line = b"Q	0	R	1	100	8M1I15M1I2M	*	*	27	GAGGATCTTGCCTGTATCGTCCGGGAG	IIIIIIIIIIIIIIIIIIIIIIIIIII"
        MIN_5 = (MATCH[0]|INS_5[0]).to_bytes()
        MIN_3 = (MATCH[0]|INS_3[0]).to_bytes()
        MIN_B = (MIN_5[0]|MIN_3[0]).to_bytes()
        expect = MATCH*6 + MIN_5 + MIN_B + MIN_3 + MATCH*11 + MIN_5 + MIN_B*2 + MIN_3 + MATCH
        muts = vectorize_read(ref, first, last, SamRead(line))
        self.assertTrue(muts == expect)


class TestVectorizeRealData(TestCase):
    """
    Test on examples of reads from real DMS-MaPseq experiments.
    """
    def test_real_1(self):
        ref = b"TATGGCTGTAGTTGTGATCAACTCCGCGAACCCATGCTTCAGTCAGCTGATGCACAATCGTTTTTAAACGGGTTTGCGGTGTAAGTGCAGCCCGTCTTAC"
        first, last = 13400, 13499
        line = b"FS10000136:77:BPG61616-2310:1:1101:1170:1610	163	SARS2	13342	44	52=1D63=1X34=	=	13474	283	CCCTGTGGGTTTTACACTTAAAAACACAGTCTGTACCGTCTGCGGTATGTGGAAGGTTATGGCTGTAGTTGTGATCAACTCCGCGAACCCATGCTTCAGTCAGCTGATGCACAATTGTTTTTAAACGGGTTTGCGGTGTAAGTGCAGCCC	FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF	AS:i:147	XN:i:0	XM:i:1	XO:i:1	XG:i:1	NM:i:2	MD:Z:52^A63C34	YS:i:149	YT:Z:CP"
        expect = MATCH*58 + SUB_T + MATCH*34 + BLANK*7
        muts = vectorize_read(ref, first, last, SamRead(line))
        self.assertTrue(muts == expect)
    
    def test_real_2(self):
        ref = b"GTCGCTTCCAAGAAAAGGACGAAGATGACAATTTAATTGATTCTTACTTTGTAGTTAAGAGACACACTTTCTCTAACTACCAACATGAAGAAACAATTTATAATTTACTTAAGGATTGTCCAGCTGTTGCTAAACATGACTTCTTTAAGTTTAGAATAGACGGTGACATGGTACCACATATATCACGTCAACGTCTTACT"
        first, last = 13601, 13800
        line = b"FS10000136:77:BPG61616-2310:1:1101:14070:2080	83	SARS2	13621	44	44=1X37=1X33=1X3=1D1X4=1I2X1=1I2=2I1X4=1X1=1I4=3S	=	13469	-297	GAAGATGACAATTTAATTGATTCTTACTTTGTAGTTAAGAGACAGACTTTCTCTAACTACCAACATGAAGAAACAATTTATATTTTACTTAAGGATTGTCCAGCTGTTGCTAAACACGACACTTTGTCGCTTCCAAGAAAAGGACGAAG	FFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF::FFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF	AS:i:119	XN:i:0	XM:i:8	XO:i:5	XG:i:6	NM:i:14	MD:Z:44C37A33T3^T0T4A0A3T4T5	YS:i:146	YT:Z:CP"
        muts = vectorize_read(ref, first, last, SamRead(line))
        # FIXME: add assertion

class TestVectorizePair(TestCase):
    def test_valid_no_muts(self):
        ref = b"ACGT"
        first, last = 1, 4
        line1 = b"Q	0	R	1	100	4=	*	*	4	ACGT	IIII"
        line2 = line1
        expect = MATCH * 4
        muts = vectorize_pair(ref, first, last, SamRead(line1), SamRead(line2))
        self.assertTrue(muts == expect)


if __name__ == "__main__":
    unittest.main()
