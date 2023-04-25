from collections import defaultdict, namedtuple
from functools import cached_property
from logging import getLogger
from pathlib import Path
import re
from typing import Iterable, Sequence

import numpy as np
import pandas as pd

from ..util.seq import DNA, BASES, A_INT, C_INT

logger = getLogger(__name__)

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


def filter_gu(seq: bytes, positions: Iterable[int]):
    """ Return each position whose base is A or C. """
    return np.asarray([pos for pos in positions
                       if seq[pos - 1] in (A_INT, C_INT)],
                      dtype=int)


def filter_polya(seq: bytes, positions: Iterable[int], max_polya: int):
    """ Discard poly(A) segments longer than the maximum length. """
    if max_polya >= 0:
        # Generate a pattern that matches poly(A) sequences longer than
        # the maximum length.
        polya_pattern = b"%c{%i,}" % (A_INT, max_polya + 1)
        # Get the 1-indexed position of every A in a poly(A) sequnce
        # longer than the maximum length.
        polya_pos = {pos for polya in re.finditer(polya_pattern, seq)
                     for pos in range(polya.start() + 1, polya.end() + 1)}
        # Remove the positions of all such poly(A) sequences.
        filtered = [pos for pos in positions if pos not in polya_pos]
    else:
        # Keep all positions.
        filtered = positions if isinstance(positions, list) else list(positions)
    # Cast to integer array.
    return np.asarray(filtered, dtype=int)


class Section(object):
    """
    Represent a section of a reference sequence between two coordinates.

    Attributes
    ----------
    ref: str
        Name of the reference sequence
    end5: int (1 ≤ end5 ≤ end3)
        Coordinate of the reference sequence at which the section's
        5' end is located (1-indexed)
    end3: int (end5 ≤ end3 ≤ len(ref_seq))
        Coordinate of the reference sequence at which the section's
        3' end is located (1-indexed; end3 itself is included)
    seq: DNA
        Sequence of the section between end5 and end3 (inclusive)
    """

    def __init__(self, /, ref_seq: DNA, ref: str, end5: int, end3: int):
        """
        Parameters
        ----------
        ref_seq: DNA
            Sequence of the entire reference (not just the section)
        ref: str
            Name of the reference sequence
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
    def name(self):
        return "full" if self.full else self.range

    @property
    def ref_name(self):
        return f"{self.ref}:{self.name}"

    def __str__(self):
        return f"Section {self.ref_name}"


class SectionFinder(Section):
    """
    The 5' and 3' ends of a section can be given explicitly as integers, but if
    the sample is of an amplicon (i.e. generated by RT-PCR using site-specific
    primers), then it is often more convenient to enter the sequences of the
    PCR primers and have the software determine the coordinates. SectionFinder
    accepts 5' and 3' coordinates given as integers or primers, validates them,
    and stores the coordinates as integers. This process works as follows:

    end5
    - If given as a parameter, this value is used.
    - Else if fwd is given, first is based on the location of its 3' end.
      - Else first is set to 1.
        - last
          - If last is given as an argument, last is set to this value.
          - Else if rev is given, last is based on the location of the 5' end
            of its reverse complement.
          - Else last is set to the length of ref_seq.

    """

    def __init__(self, /, *, ref_seq: DNA | None, ref: str, primer_gap: int,
                 end5: int | None = None, end3: int | None = None,
                 fwd: DNA | None = None, rev: DNA | None = None):
        """
        Parameters
        ----------
        ref_seq: DNA
            see superclass
        ref: DNA
            see superclass
        primer_gap: int
            (For coordinates specified by fwd/rev only) Number of
            positions 3' of the forward primer and 5' of the reverse
            primer to exclude from the section. Coordinates within 1 - 2
            nucleotides of each primer may contain DMS reactivity
            artifacts. If primer_gap = 0, then end5 and end3 are set,
            respectively, to the coordinates immediately adjacent to
            (i.e. 1 nucleotide 3' and 5' of) the 3' end of the forward
            and reverse primers.
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
        """
        if primer_gap < 0:
            logger.warning("Primer gap must be ≥ 0: setting to 0")
            primer_gap = 0
        if ref_seq is None:
            raise ValueError(f"No sequence for reference named '{ref}'. Check "
                             "that you gave the right reference sequence file "
                             "and spelled the name of the reference correctly.")
        offset = primer_gap + 1
        if end5 is None:
            # If first is to be determined from the fwd primer sequence,
            # the primer is aligned to the reference, and first is set to
            # the position (primer_gap + 1) downstream of its 3' coordinate.
            end5 = (1 if fwd is None
                    else self.locate(ref_seq, fwd).pos3 + offset)
        if end3 is None:
            # If last is to be determined from the rev primer sequence,
            # the reverse complement of the primer is aligned to the reference,
            # and last is set to the position (primer_gap + 1) upstream of its
            # 5' coordinate.
            end3 = (len(ref_seq) if rev is None
                    else self.locate(ref_seq, rev.rc).pos5 - offset)
        super().__init__(ref_seq=ref_seq, ref=ref, end5=end5, end3=end3)

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


def get_sections(ref: str, ref_seq: DNA, *,
                 coords: Iterable[tuple[int, int]],
                 primers: Iterable[tuple[DNA, DNA]],
                 primer_gap: int):
    """ Return all the sections corresponding to the given coordinates
    and/or primers in the given reference sequences. """
    sections: dict[tuple[int, int], Section] = dict()
    # Create a section for every reference and pair of coordinates.
    for end5, end3 in coords:
        try:
            section = SectionFinder(ref=ref, ref_seq=ref_seq,
                                    end5=end5, end3=end3,
                                    primer_gap=primer_gap)
        except Exception as error:
            logger.error(f"Failed to add section {ref}:{end5}-{end3}: {error}")
        else:
            if section.coord in sections:
                logger.warning(f"Skipped duplicate section: {section}")
            else:
                sections[section.coord] = section
    # Create a section for every reference and pair of primers.
    for fwd, rev in primers:
        try:
            section = SectionFinder(ref=ref, ref_seq=ref_seq,
                                    fwd=fwd, rev=rev,
                                    primer_gap=primer_gap)
        except Exception as error:
            logger.error(f"Failed to add section {ref}:{fwd}-{rev.rc}: {error}")
        else:
            if section.coord in sections:
                logger.warning(f"Skipped duplicate section: {section}")
            else:
                sections[section.coord] = section
    if not sections:
        # Create a section spanning the entire reference if no sections
        # were already defined.
        section = SectionFinder(ref=ref, ref_seq=ref_seq, primer_gap=primer_gap)
        sections[section.coord] = section
    return list(sections.values())


def get_coords_by_ref(coords: Iterable[tuple[str, int | DNA, int | DNA]]):
    ref_coords: dict[str, set[tuple[int | DNA, int | DNA]]] = defaultdict(set)
    for ref, end5, end3 in coords:
        coord = end5, end3
        if coord in ref_coords[ref]:
            logger.warning(f"Skipping duplicate coordinates: {coord}")
        else:
            ref_coords[ref].add(coord)
    return ref_coords


def sections_to_positions(sections: Iterable[Section]) -> list[int]:
    """ Return all the positions within the given (possibly overlapping)
    sections, in ascending order without duplicates. """
    return sorted(set(pos for sect in sections for pos in sect.positions))


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
                                  primer_gap=primer_gap)
            for fwd, rev in ref_primers.get(ref, list()):
                self._add_section(ref=ref, ref_seq=seq,
                                  fwd=fwd, rev=rev,
                                  primer_gap=primer_gap)
            if not self._sections[ref]:
                # If no sections were given for the reference, then add
                # one section that spans the entire reference.
                self._add_section(ref=ref, ref_seq=seq,
                                  primer_gap=primer_gap)
        if extra := (set(ref_coords) | set(ref_primers)) - set(self._sections):
            logger.warning(f"No sequences given for references: {extra}")

    def _add_section(self, **kwargs):
        """ Create a section and add it to the object. """
        try:
            section = SectionFinder(**kwargs)
        except Exception as error:
            logger.error(f"Failed to create section {kwargs}: {error}")
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
