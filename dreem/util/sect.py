from collections import Counter, defaultdict, namedtuple
from functools import cached_property
from logging import getLogger
from pathlib import Path
import re
from typing import Iterable, Sequence

import numpy as np
import pandas as pd

from .library import FIELD_REF, FIELD_START, FIELD_END, FIELD_SECT
from .seq import DNA, BASES, A_INT, C_INT

logger = getLogger(__name__)


FULL = "full"

SectionTuple = namedtuple("PrimerTuple", ["pos5", "pos3"])


def iter_ref_coords_names(library_file: Path):
    # Read every row of the library file data.
    library = pd.read_csv(library_file)
    for i, section in enumerate(zip(library[FIELD_REF], library[FIELD_START],
                                    library[FIELD_END], library[FIELD_SECT])):
        try:
            # The reference name, start coordinate, and end coordinate
            # must all have values.
            if any(map(pd.isnull, section[:-1])):
                raise ValueError(
                    f"Missing {FIELD_REF}, {FIELD_START}, and/or {FIELD_END}")
            ref = str(section[0])
            end5 = int(section[1])
            end3 = int(section[2])
            # The section name may be left blank.
            sect = "" if pd.isnull(section[3]) else str(section[3])
            # Yield the reference and coordinates as the key, and the
            # name of the section as the value.
            yield (ref, end5, end3), sect
        except Exception as error:
            logger.error(f"Failed to read section from line {i} of "
                         f"{library_file}: {error}")


def map_ref_coords_to_names(library_file: Path):
    """ Return a dictionary that maps tuples of (ref, end5, end3) to the
    names of their corresponding sections. """
    ref_coord_name_map = dict()
    for ref_coord, name in iter_ref_coords_names(library_file):
        # Check if the (ref, end5, end3) tuple has been seen already.
        if (seen := ref_coord_name_map.get(ref_coord)) is None:
            # The tuple has not been seen already: add it.
            ref_coord_name_map[ref_coord] = name
            logger.debug(
                f"Found section {ref_coord} = '{name}' in {library_file}")
        elif seen == name:
            logger.warning(f"Got duplicate section: {ref_coord} = '{name}'")
        else:
            logger.error(f"Section {ref_coord} = '{seen}' was redefined with a "
                         f"different name: '{name}' -- using name '{seen}'")
    return ref_coord_name_map


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
        raise ValueError(f"Sequence '{seq.decode()}' has {len(seq)} bases, "
                         f"but got {len(pos_list)} positions")
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
    >>> filter_polya(0, seq0)
    (b'GATCAAATCAAG', [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11])
    >>> filter_polya(1, seq0)
    (b'GTCTCG', [0, 2, 3, 7, 8, 11])
    >>> filter_polya(1, seq0, positions=range(10, 22))
    (b'GTCTCG', [10, 12, 13, 17, 18, 21])
    >>> filter_polya(3, seq0)
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
                 end5: int | None = None, end3: int | None = None,
                 name: str | None = None):
        """
        Parameters
        ----------
        ref: str
            Name of the reference sequence
        ref_seq: DNA
            Sequence of the entire reference (not just the section)
        end5: int | None = None
            Coordinate of the reference sequence at which the 5' end of
            the section is located. If positive, number the coordinates
            of the reference sequence 1, 2, ... starting at the 5' end
            (i.e. 1-based indexing). If negative, number the coordinates
            of the reference sequence -1, -2, ... starting at the 3' end
            (i.e. 1-based indexing from the other side), then convert to
            the corresponding (positive) 1-based index from the 5' end
        end3: int | None = None
            Coordinate of the reference sequence at which the section's
            3' end is located. Follows the same coordinate numbering
            convention as end5
        name: str | None = None
            Name of the section. If blank, defaults to ```self.range```.
        """
        self.ref = ref
        if end5 is None and end3 is None:
            # If both coordinates are omitted, use the full reference
            # and set the name to full.
            end5 = 1
            end3 = len(ref_seq)
            name = FULL
        if end5 < 0:
            # Compute the corresponding positive coordinate.
            end5 += len(ref_seq) + 1
        if end3 < 0:
            # Compute the corresponding positive coordinate.
            end3 += len(ref_seq) + 1
        self.end5: int = end5
        self.end3: int = end3
        if not 1 <= end5 <= end3 <= len(ref_seq):
            raise ValueError("Must have 1 ≤ end5 ≤ end3 ≤ len(ref_seq), "
                             f"but got end5 = {end5}, end3 = {end3}, and "
                             f"len(ref_seq) = {len(ref_seq)}")
        self.seq: DNA = ref_seq[end5 - 1: end3]
        self.full = end5 == 1 and end3 == len(ref_seq)
        if name is None:
            self.name = self.range
        elif isinstance(name, str):
            self.name = name if name else self.range
        else:
            raise TypeError(f"Parameter 'name' must be 'str', not {type(str)}")

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

    def to_dict(self):
        return {"ref": self.ref, "seq": self.seq, "sect": self.name,
                "end5": self.end5, "end3": self.end3}

    def __str__(self):
        if self.name == self.range:
            return f"Section {self.ref_name}"
        else:
            return f"Section {self.ref_name} ({self.end5}-{self.end3})"


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
        # Get the names of the sections from the library, if any.
        section_names = dict()
        if library is not None:
            try:
                section_names = map_ref_coords_to_names(library)
                # Combines the coordinates from the library and from the
                # coord parameter.
                coords = list(coords) + list(section_names)
            except Exception as error:
                logger.error(f"Failed to add coordinates from library: {error}")

        # Group coordinates and primers by reference.
        ref_coords = get_coords_by_ref(coords)
        ref_primers = get_coords_by_ref(primers)

        # For each reference, generate sections from the coordinates.
        self._sections: dict[str, dict[tuple[int, int], Section]] = dict()
        for ref, seq in ref_seqs:
            self._sections[ref] = dict()
            for end5, end3 in ref_coords[ref]:
                # Add a section for each pair of 5' and 3' coordinates.
                self._add_section(ref=ref, ref_seq=seq, end5=end5, end3=end3,
                                  primer_gap=primer_gap,
                                  name=section_names.get((ref, end5, end3)))
            for fwd, rev in ref_primers[ref]:
                # Add a section for each pair of fwd and rev primers.
                self._add_section(ref=ref, ref_seq=seq, fwd=fwd, rev=rev,
                                  primer_gap=primer_gap)
            if not self._sections[ref]:
                # If no sections were given for the reference, then add
                # a section named 'full' that spans the full reference.
                self._add_section(ref=ref, ref_seq=seq,
                                  primer_gap=primer_gap)
        if extra := (set(ref_coords) | set(ref_primers)) - set(self._sections):
            logger.warning(f"No sequences given for references: {extra}")

    def _add_section(self, **kwargs):
        """ Create a section and add it to the object. """
        try:
            section = SectionFinder(**kwargs)
        except Exception as error:
            logger.error(f"Failed to create section with {kwargs}: {error}")
        else:
            # Check if the section was seen already.
            if (seen := self._sections[section.ref].get(section.coord)) is None:
                # The section was not seen already: add it.
                self._sections[section.ref][section.coord] = section
            elif seen.name == section.name:
                # The section was seen already with the same name.
                logger.warning(f"Got duplicate section: {section}")
            else:
                # The section was seen already with a different name.
                logger.error(f"Section {seen} was redefined with as {section}. "
                             f"Using the first encountered: {seen}")

    def list(self, ref: str):
        """ Return a list of the sections for a given reference. """
        return list(self._sections[ref].values())

    @property
    def count(self):
        """ Total number of sections. """
        return sum(map(len, self._sections.values()))
