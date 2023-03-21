from collections import namedtuple
import os
from typing import Any, Iterable

import pandas as pd
import pytest

import dreem
from ...vector.profile import get_min_qual
from ...vector.vector import SamRead, vectorize_read


DREEM_DIR = os.path.dirname(os.path.abspath(dreem.__file__))
DATA_FILE = os.path.join(DREEM_DIR, "test_files", "vectorize_read_tests.csv")

SAM_FIELDS = ("QNAME", "FLAG", "RNAME", "POS", "MAPQ", "CIGAR", "RNEXT",
              "PNEXT", "TLEN", "SEQ", "QUAL")
SAM_LINE_TEMPLATE = "\t".join("{" + field + "}" for field in SAM_FIELDS)

HEX_PREFIX = "0x"
HEX_LENGTH = 2
HEX_BASE = 16


test_result = namedtuple("test_result", ("name", "desc", "expect", "result"))


def hex_to_bytes(hex_str: str):
    """ Convert a hexadecimal string to an array of bytes. """
    if not hex_str.startswith(HEX_PREFIX):
        raise ValueError(f"Not a hexadecimal string: '{hex_str}'")
    # Compute the number of bytes in the hexadecimal string.
    n_bytes, rem = divmod(len(hex_str) - len(HEX_PREFIX), HEX_LENGTH)
    if rem:
        raise ValueError(f"Length of hexadecimal string '{hex_str}' "
                         f"is not a multiple of {HEX_LENGTH}")
    # Convert each set of hexadecimal characters to an integer, and then
    # assemble them into a bytes object (implicitly converting to bytes)
    return bytes(int(hex_str[i: i + HEX_LENGTH], HEX_BASE)
                 for i in range(len(HEX_PREFIX),
                                len(HEX_PREFIX) + n_bytes * HEX_LENGTH,
                                HEX_LENGTH))


def get_sam_read(df: pd.DataFrame, index: Any):
    """ Return a SamRead instance from a row in a DataFrame. """
    # Collect SAM fields from dataframe.
    sam_fields = {field: df.loc[index, field] for field in SAM_FIELDS}
    # Create the line in SAM format.
    sam_line = SAM_LINE_TEMPLATE.format(**sam_fields)
    # Return a SamRead from the line.
    return SamRead(sam_line.encode())


def run_row(df: pd.DataFrame, index: Any):
    """ Test one row from the dataframe. """
    # Create a SamRead from the row.
    read = get_sam_read(df, index)
    # Gather test parameters from the row.
    name = str(df.loc[index, "Test"])
    desc = str(df.loc[index, "Description"])
    rseq = str(df.loc[index, "Section"]).encode()
    end5 = int(df.loc[index, "End5"])
    end3 = int(df.loc[index, "End3"])
    penc = int(df.loc[index, "PhredEnc"])
    pmin = int(df.loc[index, "MinPhred"])
    ambid = bool(df.loc[index, "Ambid"])
    expect = hex_to_bytes(df.loc[index, "Expected"])
    qmin = get_min_qual(pmin, penc)
    # Vectorize the read.
    result = vectorize_read(read,
                            sect_seq=rseq,
                            section_end5=end5,
                            min_qual=qmin,
                            ambid=ambid)
    # Return the results.
    return test_result(name=name, desc=desc, expect=expect, result=result)


def compile_results(results: Iterable[test_result]):
    df = pd.DataFrame.from_records(results, columns=test_result._fields)
    df.set_index("name", drop=True)
    return df


def get_pass_fail(results: pd.DataFrame):
    matched = results["result"] == results["expect"]
    passed = list(results.loc[matched, "name"])
    failed = list(results.loc[~matched, "name"])
    return passed, failed


@pytest.mark.skip(reason="Function must be called via test_run")
def run_all_rows(df: pd.DataFrame):
    results = compile_results(run_row(df, index) for index in df.index)
    return get_pass_fail(results)


def run_csv_file(csv_file):
    passed, failed = run_all_rows(pd.read_csv(csv_file))
    print(f"Passed (n = {len(passed)}): {passed}")
    print(f"Failed (n = {len(failed)}): {failed}")


@pytest.mark.skipif(not os.path.isfile(DATA_FILE),
                    reason=f"Data file not found: {DATA_FILE}")
def test_run():
    run_csv_file(DATA_FILE)
