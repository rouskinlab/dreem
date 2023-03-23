import numpy as np
import pandas as pd

from dreem.vector.profile import VectorReader


N_BITS = 8
BYTE_LIMIT = 2**N_BITS
POWERS = np.power(2, np.arange(N_BITS, dtype=int))


def get_all_bytes():
    return range(BYTE_LIMIT)


def get_vector_df(vector: int):
    assert 0 <= vector < BYTE_LIMIT
    return pd.DataFrame(np.full(shape=(1, 1),
                                fill_value=vector,
                                dtype=np.uint8))


def byte_to_bits(byte: int):
    """ Return an 8-element boolean array of the bits in a byte. """
    assert 0 <= byte < BYTE_LIMIT
    return np.asarray(np.bitwise_and(byte, POWERS), dtype=bool)


def vector_is_subset(vector: int, query: int):
    """ Return whether vector is a bitwise subset of query. """
    if vector == 0:
        # Return False if vector is blank.
        return False
    vector_bits = byte_to_bits(vector)
    query_bits = byte_to_bits(query)
    # Vector is a subset if there are not any bits that are
    # 1 in vector and 0 in query.
    return not np.any(np.logical_and(vector_bits,
                                     np.logical_not(query_bits)
                                     ))


def vector_is_superset(vector: int, query: int):
    """ Return whether vector is a bitwise superset of query. """
    if vector == 0:
        # Return False if vector is blank.
        return False
    vector_bits = byte_to_bits(vector)
    query_bits = byte_to_bits(query)
    # Vector is a superset if there are not any bits that are
    # 1 in query and 0 in vector.
    return not np.any(np.logical_and(query_bits,
                                     np.logical_not(vector_bits)
                                     ))


def vector_is_match(vector: int, query: int):
    """ Return whether vector matches query. """
    # Return False if vector is blank.
    return (vector == query) and (vector != 0)


def check_vector_query(vector: int, query: int):
    """ Check if the result of VectorReader._query_vectors equals the
    true result. """
    is_match = vector_is_match(vector, query)
    is_subset = vector_is_subset(vector, query)
    is_superset = vector_is_superset(vector, query)
    expect = np.array([is_match,
                       is_subset,
                       is_superset,
                       is_subset or is_superset])
    vector_df = get_vector_df(vector)
    res_match = np.all(VectorReader._query_vectors(vector_df, query))
    res_subset = np.all(VectorReader._query_vectors(vector_df, query,
                                                    subsets=True))
    res_superset = np.all(VectorReader._query_vectors(vector_df, query,
                                                      supersets=True))
    res_both = np.all(VectorReader._query_vectors(vector_df, query,
                                                  subsets=True,
                                                  supersets=True))
    result = np.array([res_match, res_subset, res_superset, res_both])
    return expect, result


def test_run(verbose: bool = False):
    succeeded = list()
    failed = list()
    for query in get_all_bytes():
        print(query)
        for vector in get_all_bytes():
            case = query, vector
            expect, result = check_vector_query(vector, query)
            success = np.all(expect == result)
            if verbose:
                print(f"{case}\tExp: {expect}\tRes: {result}\tPass: {success}")
            if success:
                succeeded.append(case)
            else:
                failed.append(case)
    print("Succeeded:", len(succeeded))
    print("Failed:", len(failed))


if __name__ == "__main__":
    test_run()
