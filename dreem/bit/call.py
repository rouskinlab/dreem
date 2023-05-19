import numpy as np
import pandas as pd

from ..util.sect import cols_to_seq_pos
from ..util.seq import (MATCH, DELET, INS_5, INS_3,
                        SUB_A, SUB_C, SUB_G, SUB_T,
                        A_INT, C_INT, G_INT, T_INT)


class BitCaller(object):
    def __init__(self, *,
                 subset: bool = False,
                 aa: bool = False,
                 ac: bool = False,
                 ag: bool = False,
                 at: bool = False,
                 ad: bool = False,
                 ai: bool = False,
                 ca: bool = False,
                 cc: bool = False,
                 cg: bool = False,
                 ct: bool = False,
                 cd: bool = False,
                 ci: bool = False,
                 ga: bool = False,
                 gc: bool = False,
                 gg: bool = False,
                 gt: bool = False,
                 gd: bool = False,
                 gi: bool = False,
                 ta: bool = False,
                 tc: bool = False,
                 tg: bool = False,
                 tt: bool = False,
                 td: bool = False,
                 ti: bool = False):
        self._sub = subset
        self._queries = {A_INT: self._make_query(aa, False, ac, ag, at, ad, ai),
                         C_INT: self._make_query(cc, ca, False, cg, ct, cd, ci),
                         G_INT: self._make_query(gg, ga, gc, False, gt, gd, gi),
                         T_INT: self._make_query(tt, ta, tc, tg, False, td, ti)}

    @staticmethod
    def _make_query(nn: bool, _a: bool, _c: bool, _g: bool, _t: bool,
                    _d: bool, _i: bool):
        query = (MATCH | INS_5) if nn else 0
        query |= SUB_A if _a else 0
        query |= SUB_C if _c else 0
        query |= SUB_G if _g else 0
        query |= SUB_T if _t else 0
        query |= DELET if _d else 0
        query |= INS_3 if _i else 0
        return query

    def _call_df(self, mut_vectors: pd.DataFrame, query: int) -> pd.DataFrame:
        return np.equal(query, (np.bitwise_or(mut_vectors, query) if self._sub
                                else mut_vectors))

    def _call_base(self,
                   bits: pd.DataFrame,
                   mut_vectors: pd.DataFrame,
                   seq: np.ndarray,
                   base: int):
        # Determine the positions in the sequence that have the base.
        pos = np.equal(seq, base)
        if not np.any(pos):
            # If no bases match, there is nothing to do.
            return
        # Set the bits of the positions that have the base.
        bits.loc[:, pos] = self._call_df(mut_vectors.loc[:, pos],
                                         self._queries[base])

    def call(self, mut_vectors: pd.DataFrame):
        """
        Query the given mutation vectors. Return a boolean DataFrame
        indicating which elements of mutation vectors matched the query.

        Parameters
        ----------
        mut_vectors: DataFrame
            Mutation vectors

        Returns
        -------
        DataFrame
            Boolean DataFrame of the same shape as mut_vectors wherein
            element i,j is True if element i,j of mut_vectors matches
            the query, otherwise False.
        """
        # Get the reference sequence of the mutation vectors as an array
        # of 8-bit unsigned integers.
        seq = np.frombuffer(cols_to_seq_pos(mut_vectors.columns.to_list())[0],
                            dtype=np.uint8)
        # Initialize all bits to False (i.e. not matching query).
        bits = pd.DataFrame(np.zeros_like(mut_vectors),
                            index=mut_vectors.index,
                            columns=mut_vectors.columns)
        # Call the bits for each type of base.
        for base in self._queries:
            print("Calling base", chr(base))
            self._call_base(bits, mut_vectors, seq, base)
        return bits

    @classmethod
    def mut(cls, count_del: bool, count_ins: bool):
        """ Return a BitCaller that marks mutations as True. """
        return cls(ac=True, ag=True, at=True, ad=count_del, ai=count_ins,
                   ca=True, cg=True, ct=True, cd=count_del, ci=count_ins,
                   ga=True, gc=True, gt=True, gd=count_del, gi=count_ins,
                   ta=True, tc=True, tg=True, td=count_del, ti=count_ins,
                   subset=True)

    @classmethod
    def ref(cls):
        """ Return a BitCaller that marks matches as True. """
        return cls(aa=True, cc=True, gg=True, tt=True, subset=True)

    @classmethod
    def cov(cls):
        """ Return a BitCaller that marks covered positions as True. """
        return cls(aa=True, ac=True, ag=True, at=True, ad=True, ai=True,
                   ca=True, cc=True, cg=True, ct=True, cd=True, ci=True,
                   ga=True, gc=True, gg=True, gt=True, gd=True, gi=True,
                   ta=True, tc=True, tg=True, tt=True, td=True, ti=True,
                   subset=True)
