from hashlib import md5
from typing import Any


def digest_file(file_path: Any) -> str:
    """
    Compute the checksum of a file.

    Parameters
    ----------
    file_path: Any (path-like)
        Path of the file on which to compute the checksum. Can be
        any type that the open() function recognizes as a path.

    Returns
    -------
    str
        Checksum of the file (in hexadecimal)
    """
    with open(file_path, "rb") as f:
        digest = md5(f.read()).hexdigest()
    return digest
