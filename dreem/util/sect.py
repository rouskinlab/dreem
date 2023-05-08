from collections import Counter, defaultdict, namedtuple
from functools import cached_property
from logging import getLogger
from pathlib import Path
import re
from typing import Iterable, Sequence

import numpy as np
import pandas as pd

from ..util.seq import DNA, BASES, A_INT, C_INT

logger = getLogger(__name__)


FULL = "full"

SectionTuple = namedtuple("PrimerTuple", ["pos5", "pos3"])


def add_coords_from_library(library_path: str,
                            initial_coords: tuple[tuple[str, int, int]]):
    library_coords = list()
    try:
        library = pd.read_csv(library_path)
        for ref, end5, end3 in zip(library["reference"],
                                   library["section_start"],
                                   library["section_end"], strict=True):
            try:
                coord = (str(Path(ref).with_suffix("")),
                         int(end5),
                         int(end3))
            except ValueError as error:
                logger.error(f"Failed to add coordinates {ref, end5, end3} "
                             f"with the following error: {error}")
            else:
                if coord in initial_coords or coord in library_coords:
                    logger.warning(f"Skipping duplicate coordinates: {coord}")
                else:
                    library_coords.append(coord)
    except (FileNotFoundError, KeyError, ValueError) as error:
        logger.critical(f"Failed to add coordinates from {library_path} "
                        f"with the following error: {error}")
    return initial_coords + tuple(library_coords)


def encode_primer(ref: str, fwd: str, rev: str):
    """ Convert a pair of primers from strings to DNA objects. """
    return ref, DNA(fwd.encode()), DNA(rev.encode())


def encode_primers(primers: Iterable[tuple[str, str, str]]):
    """ Convert pairs of primers from strings to DNA objects. """
    enc_primers = dict()
    for primer in primers:
        if primer in enc_primers:
            logger.warning(f"Skipping duplicate primer: {primer}")
        else:
            try:
                enc_primers[primer] = encode_primer(*primer)
            except Exception as error:
                logger.error(f"Failed to encode primer {primer}: {error}")
                enc_primers[primer] = None
    # Return a list of all encoded primers with None values removed.
    return list(filter(None, enc_primers.values()))


def seq_pos_to_cols(seq: bytes, positions: Iterable[int]):
    """ Convert sequence and positions to column names. Each column
    name is a string of the base followed by the position. """
    # Use chr(base) instead of base.decode() because base is an int.
    # Subtract 1 from pos when indexing seq because pos is 1-indexed.
    return [f"{chr(seq[pos - 1])}{pos}" for pos in positions]


def cols_to_seq_pos(columns: list[str]):
    """ Convert column names to sequence and positions. Each column
    name is a string of the base followed by the position. """
    # Regex pattern "^([ACGT])([0-9]+)$" finds the base and position
    # in each column name.
    pattern = re.compile(f"^([{BASES.decode()}])([0-9]+)$")
    # Match each column name using the pattern.
    matches = list(map(pattern.match, columns))
    try:
        # Obtain the two groups (base and position) from each match
        # and unzip them into two tuples.
        bases, pos_strs = zip(*map(re.Match.groups, matches))
    except TypeError:
        # TypeError is raised if any match is None, which happens if
        # a column fails to match the pattern.
        invalid_cols = [col for col, match in zip(columns, matches,
                                                  strict=True)
                        if match is None]
        raise ValueError(f"Invalid columns: {invalid_cols}")
    # Join the tuple of bases into a DNA sequence.
    seq = DNA("".join(bases).encode())
    # Cast the tuple of strings of positions into an integer array.
    positions = np.asarray(list(map(int, pos_strs)), dtype=int)
    return seq, positions


def positions_from_seq(seq: bytes,
                       positions: Iterable[int] | None) -> list[int]:
    """ Return positions for a sequence, optionally from given
    positions. """
    if positions is None:
        # If no positions were given, return [0, 1, ... , len(seq) - 1].
        return list(range(len(seq)))
    # Check for duplicate positions.
    pos_counts = Counter(positions)
    if pos_dups := [pos for pos, count in pos_counts.items() if count > 1]:
        raise ValueError(f"Got duplicate positions: {pos_dups}")
    # List the positions (in their original order) from the counts.
    pos_list = list(pos_counts.keys())
    # Verify that there are as many positions as bases in seq.
    if len(pos_list) != len(seq):
        raise ValueError(f"Sequence has {len(seq)} bases, but got "
                         f"{len(pos_list)} positions")
    return pos_list


def filter_gu(seq: bytes, positions: Iterable[int] | None = None):
    """ Return each position whose base is A or C.

    Parameters
    ----------
    seq: bytes
        Sequence from which to remove bases other than A and C.
    positions: Iterable[int] | None = None
        Numeric positions of the sequence. If given, must be an iterable
        of unique integers of the same length as ```seq```. If omitted,
        defaults to ```range(len(seq))```.

    Returns
    -------
    tuple[bytes, list[int]]
        List of the integer positions from ```positions``` that remain
        after removing bases other than A and C.
    """
    # Find the positions excluding those adenines.
    idx_out, pos_out = zip(*[(idx, pos) for idx, pos
                             in enumerate(positions_from_seq(seq, positions))
                             if seq[idx] in (A_INT, C_INT)])
    # Also generate the sequence without the Gs and Us.
    seq_out = bytes(map(seq.__getitem__, idx_out))
    return seq_out, list(pos_out)


def filter_polya(exclude_polya: int,
                 seq: bytes,
                 positions: Iterable[int] | None = None):
    """ Discard poly(A) stretches with length ≥ ```exclude_polya```.

    Parameters
    ----------
    exclude_polya: int
        Remove all positions in ```seq``` that are part of a stretch of
        consecutive A bases of length ≥ ```exclude_polya```. If 0, then do
        not remove any positions. Must be ≥ 0.
    seq: bytes
        Sequence from which to remove consecutive A bases.
    positions: Iterable[int] | None = None
        Numeric positions of the sequence. If given, must be an iterable
        of unique integers of the same length as ```seq```. If omitted,
        defaults to ```range(len(seq))```.

    Returns
    -------
    tuple[bytes, list[int]]
        - [0] Sequence with poly(A) sequences removed
        - [1] List of remaining positions.

    Examples
    --------
    >>> seq0 = b"GATCAAATCAAG"
    >>> filter_polya(seq0, 0)
    (b'GATCAAATCAAG', [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11])
    >>> filter_polya(seq0, 1)
    (b'GTCTCG', [0, 2, 3, 7, 8, 11])
    >>> filter_polya(seq0, 1, positions=range(10, 22))
    (b'GTCTCG', [10, 12, 13, 17, 18, 21])
    >>> filter_polya(seq0, 3)
    (b'GATCTCAAG', [0, 1, 2, 3, 7, 8, 9, 10, 11])
    """
    if exclude_polya < 0:
        raise ValueError(f"exclude_polya must be ≥ 0, but got {exclude_polya}")
    if exclude_polya == 0:
        # Do not remove any positions.
        return bytes(seq), positions_from_seq(seq, positions)
    # Generate a pattern that matches stretches of consecutive adenines
    # that are at least as long as exclude_polya.
    polya_pattern = b"%c{%i,}" % (A_INT, exclude_polya)
    # Get the 0-indexed position of every adenine in those stretches.
    polya_idxs = {pos for polya in re.finditer(polya_pattern, seq)
                  for pos in range(polya.start(), polya.end())}
    # Find the positions excluding those adenines.
    idx_out, pos_out = zip(*[(idx, pos) for idx, pos
                             in enumerate(positions_from_seq(seq, positions))
                             if idx not in polya_idxs])
    # Also generate the sequence without those adenines.
    seq_out = bytes(map(seq.__getitem__, idx_out))
    return seq_out, list(pos_out)


def filter_pos(exclude_pos: Iterable[int],
               seq: bytes,
               positions: Iterable[int] | None = None):
    """ Discard arbitrary positions given in ```exclude_pos```. """
    exclude_pos_set = set(exclude_pos)
    # Find the non-excluded positions.
    idx_out, pos_out = zip(*[(idx, pos) for idx, pos
                             in enumerate(positions_from_seq(seq, positions))
                             if pos not in exclude_pos_set])
    # Also generate the sequence without those positions.
    seq_out = bytes(map(seq.__getitem__, idx_out))
    return seq_out, list(pos_out)


class Section(object):
    """
    Represent a section of a reference sequence between two coordinates.
    """

    def __init__(self, /, ref: str, ref_seq: DNA, *,
                 name: str | None = None, end5: int, end3: int):
        """
        Parameters
        ----------
        ref: str
            Name of the reference sequence
        ref_seq: DNA
            Sequence of the entire reference (not just the section)
        name: str | None = None
            Name of the section. If falsy, defaults to ```self.range```.
        end5: int (-len(ref_seq) ≤ end5 ≤ len(ref_seq); end5 ≠ 0)
            Coordinate of the reference sequence at which the 5' end of
            the section is located. If positive, number the coordinates
            of the reference sequence 1, 2, ... starting at the 5' end
            (i.e. 1-based indexing). If negative, number the coordinates
            of the reference sequence -1, -2, ... starting at the 3' end
            (i.e. 1-based indexing from the other side), then convert to
            the corresponding (positive) 1-based index from the 5' end
        end3: int (-len(ref_seq) ≤ end3 ≤ len(ref_seq); end3 ≠ 0)
            Coordinate of the reference sequence at which the section's
            3' end is located. Follows the same coordinate numbering
            convention as end5
        """
        self.ref = ref
        if end5 < 0:
            # Compute the corresponding positive coordinate.
            end5 += len(ref_seq) + 1
        if end3 < 0:
            # Compute the corresponding positive coordinate.
            end3 += len(ref_seq) + 1
        self.end5 = end5
        self.end3 = end3
        if not 1 <= end5 <= end3 <= len(ref_seq):
            raise ValueError("Must have 1 ≤ end5 ≤ end3 ≤ len(ref_seq), "
                             f"but got end5 = {end5}, end3 = {end3}, and "
                             f"len(ref_seq) = {len(ref_seq)}")
        self.seq = ref_seq[end5 - 1: end3]
        self.full = end5 == 1 and end3 == len(ref_seq)
        self.name = str(name) if name else self.range

    @property
    def length(self):
        """ Return the length of the section of interest. """
        return self.end3 - self.end5 + 1

    @property
    def coord(self):
        """ Return the 5' and 3' coordinates as a tuple. """
        return self.end5, self.end3

    @cached_property
    def positions(self):
        """ Return all positions in the section of interest. """
        return np.arange(self.end5, self.end3 + 1, dtype=int)

    def subseq(self, positions: Sequence[int] | None):
        """ Return a subset of the sequence at the given positions, or
        the entire sequence if positions is None. """
        if positions is None:
            # Return the entire sequence if no positions are selected.
            return self.seq
        n_pos = len(positions)
        if n_pos == 0:
            raise ValueError("Positions is an empty sequence")
        pos5, pos3 = positions[0], positions[-1]
        if n_pos != pos3 - pos5 + 1:
            raise ValueError(
                "Positions must be a sequence of monotonically increasing "
                f"consecutive integers, but got {positions}")
        if pos5 < self.end5:
            raise ValueError(f"5' end ({pos5}) out of bounds for {self}")
        if pos3 > self.end3:
            raise ValueError(f"3' end ({pos3}) out of bounds for {self}")
        return self.seq[pos5 - self.end5: pos3 - self.end5 + 1]

    @property
    def ref_coord(self):
        """ Return the name of the reference and the 5' and 3' positions
        of the section of interest; for hashing and equality test_input. """
        return self.ref, self.end5, self.end3

    @cached_property
    def columns(self):
        """ Return the column names of the section. """
        return seq_pos_to_cols(self.seq, self.positions)

    @property
    def range(self):
        return f"{self.end5}-{self.end3}"

    @property
    def ref_name(self):
        return f"{self.ref}:{self.name}"

    def __str__(self):
        return f"Section {self.ref_name}"


class SectionFinder(Section):
    """
    The 5' and 3' ends of a section can be given explicitly as integers,
    but if the sample is of an amplicon (i.e. generated by RT-PCR using
    site-specific primers), then it is often more convenient to enter
    the sequences of the PCR primers and have the software determine the
    coordinates. SectionFinder accepts 5' and 3' coordinates given as
    integers or primers, validates them, and stores the coordinates as
    integers, as follows:

    end5 = end5 if end5 is given, else the 3' end of the forward primer
           + (primer_gap + 1) if fwd is given, else 1
    end3 = end3 if end3 is given, else the 5' end of the reverse primer
       - (primer_gap + 1) if rev is given, else the length of ref_seq
    """

    def __init__(self, /, ref: str, ref_seq: DNA | None, *,
                 name: str | None = None,
                 end5: int | None = None, end3: int | None = None,
                 fwd: DNA | None = None, rev: DNA | None = None,
                 primer_gap: int | None = None):
        """
        Parameters
        ----------
        ref: str
            see superclass
        ref_seq: DNA
            see superclass
        name: str | None = None
            see superclass
        end5: int | None = None
            If given, behaves as in the superclass; otherwise, ignored.
        end3: int | None = None
            If given, behaves as in the superclass; otherwise, ignored.
        fwd: DNA | None = None
            (For amplicons only) Sequence of the forward PCR primer
            that was used to generate the amplicon
        rev: DNA | None = None
            (For amplicons only) Sequence of the reverse PCR primer
            that was used to generate the amplicon (the actual sequence,
            not its reverse complement)
        primer_gap: int | None = None
            (For coordinates specified by fwd/rev only) Number of
            positions 3' of the forward primer and 5' of the reverse
            primer to exclude from the section. Coordinates within 1 - 2
            nucleotides of each primer may contain DMS reactivity
            artifacts. If primer_gap = 0, then end5 and end3 are set,
            respectively, to the coordinates immediately adjacent to
            (i.e. 1 nucleotide 3' and 5' of) the 3' end of the forward
            and reverse primers.
        """
        if primer_gap < 0:
            raise ValueError(f"primer_gap must be ≥ 0, but got {primer_gap}")
        if ref_seq is None:
            raise ValueError(f"No sequence for reference named '{ref}'. Check "
                             "that you gave the right reference sequence file "
                             "and spelled the name of the reference correctly.")
        if end5 is None:
            # No 5' end coordinate was given.
            if fwd is None:
                # No forward primer was given: default to first base.
                end5 = 1
            elif primer_gap is None:
                raise TypeError("Must give primer_gap if using primers")
            else:
                # Locate the end of the forward primer, then end the
                # section (primer_gap + 1) positions downstream.
                end5 = self.locate(ref_seq, fwd).pos3 + (primer_gap + 1)
        if end3 is None:
            # No 3' end coordinate was given.
            if rev is None:
                # No reverse primer was given: default to last base.
                end3 = len(ref_seq)
            elif primer_gap is None:
                raise TypeError("Must give primer_gap if using primers")
            else:
                # Locate the start of the reverse primer, then end the
                # section (primer_gap + 1) positions upstream.
                end3 = self.locate(ref_seq, rev.rc).pos5 - (primer_gap + 1)
        super().__init__(ref_seq=ref_seq, ref=ref,
                         end5=end5, end3=end3, name=name)

    @staticmethod
    def locate(ref_seq: DNA, primer: DNA) -> SectionTuple:
        """
        Return the 5' and 3' positions (1-indexed) of a primer within a
        reference sequence. The primer must occur exactly once in the
        reference, otherwise an error is raised.

        Parameters
        ----------
        ref_seq: DNA
            Sequence of the entire reference (not just the section of
            interest)
        primer: DNA
            Sequence of the forward PCR primer or the reverse complement
            of the reverse PCR primer

        Returns
        -------
        SectionTuple
            Named tuple of the first and last positions that the primer
            occupies in the reference sequence. Positions are 1-indexed
            and include the first and last coordinates.
        """
        matches = list(re.finditer(primer, ref_seq))
        if not matches:
            raise ValueError(f"Primer '{primer}' is not in ref '{ref_seq}'")
        if len(matches) > 1:
            raise ValueError(f"Primer '{primer}' occurs {len(matches)} times "
                             f"in ref '{ref_seq}'")
        # Add 1 to convert from 0-indexed to 1-indexed.
        pos5 = matches[0].start() + 1
        # No change is needed to convert from exclusive 0-indexed to
        # inclusive 1-indexed.
        pos3 = matches[0].end()
        return SectionTuple(pos5, pos3)


def sects_to_pos(sections: Iterable[Section]) -> list[int]:
    """ Return all the positions within the given (possibly overlapping)
    sections, in ascending order without duplicates. """
    return sorted(set(pos for sect in sections for pos in sect.positions))


def get_coords_by_ref(coords: Iterable[tuple[str, int | DNA, int | DNA]]):
    ref_coords: dict[str, set[tuple[int | DNA, int | DNA]]] = defaultdict(set)
    for ref, end5, end3 in coords:
        coord = end5, end3
        if coord in ref_coords[ref]:
            logger.warning(f"Skipping duplicate coordinates: {coord}")
        else:
            ref_coords[ref].add(coord)
    return ref_coords


class RefSections(object):
    """ A collection of sections, grouped by reference. """

    def __init__(self,
                 ref_seqs: Iterable[tuple[str, DNA]], *,
                 library: Path | None = None,
                 coords: Iterable[tuple[str, int, int]] = (),
                 primers: Iterable[tuple[str, DNA, DNA]] = (),
                 primer_gap: int):
        # Process the library.
        if library and False:  # FIXME: add library
            df_library = pd.read_csv(library)
            section_names = get_library_sections(df_library)
            if section_names:
                # Add the coordinates from the library to the coordinates
                # given as a parameter. Remove duplicates and preserve the
                # order by converting the coordinates into dictionary keys
                # first (with Counter) and then casting back to a tuple.
                coords = tuple(Counter(coords + tuple(section_names)))
        else:
            section_names = dict()

        # Group coordinates and primers by reference.
        ref_coords = get_coords_by_ref(coords)
        ref_primers = get_coords_by_ref(primers)

        # For each reference, generate sections from the coordinates.
        self._sections: dict[str, dict[tuple[int, int], Section]] = dict()
        for ref, seq in ref_seqs:
            self._sections[ref] = dict()
            for end5, end3 in ref_coords.get(ref, list()):
                self._add_section(ref=ref, ref_seq=seq,
                                  end5=end5, end3=end3,
                                  primer_gap=primer_gap,
                                  name=section_names.get((end5, end3)))
            for fwd, rev in ref_primers.get(ref, list()):
                self._add_section(ref=ref, ref_seq=seq,
                                  fwd=fwd, rev=rev,
                                  primer_gap=primer_gap)
            if not self._sections[ref]:
                # If no sections were given for the reference, then add
                # a section named 'full' that spans the full reference.
                self._add_section(ref=ref, ref_seq=seq,
                                  primer_gap=primer_gap,
                                  name=FULL)
        if extra := (set(ref_coords) | set(ref_primers)) - set(self._sections):
            logger.warning(f"No sequences given for references: {extra}")

    def _add_section(self, **kwargs):
        """ Create a section and add it to the object. """
        try:
            section = SectionFinder(**kwargs)
        except Exception as error:
            logger.error(f"Failed to create section with {kwargs}: {error}")
        else:
            if section.coord in self._sections[section.ref]:
                logger.warning(f"Duplicate section: {section}")
            else:
                self._sections[section.ref][section.coord] = section

    def list(self, ref: str):
        return list(self._sections[ref].values())

    @property
    def count(self):
        """ Total number of sections. """
        return sum(map(len, self._sections.values()))
