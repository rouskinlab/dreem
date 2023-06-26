"""
Core -- Sections Module
========================================================================
Auth: Matty

Utilities for sections of reference sequences.
"""

from collections import defaultdict, namedtuple
from functools import cached_property, reduce
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


# FIXME: Eventually replace the dual seq-pos and pos indexing schemes
# FIXME: with one indexing scheme that uses a MultiIndex with the
# FIXME: positions at level 0 and the bases at level 1.

def seq_pos_to_index(seq: bytes, positions: Sequence[int], start: int):
    """
    Convert sequence and positions to indexes, where each index is a
    string of the base followed by the position.

    Parameters
    ----------
    seq: bytes
        DNA or RNA sequence.
    positions: Sequence[int]
        Positions of the sequence from which to build the index.
    start: int
        Numerical position to assign to the first base in the sequence.

    Returns
    -------
    pd.Index
        Index of the same length as positions where each element is a
        string of {base}{position}.
    """
    return pd.Index([f"{chr(seq[pos - start])}{pos}" for pos in positions])


def index_to_seq_pos(index: pd.Index):
    """ Convert index names to sequence and positions. Each index name
    is a string of the base followed by the position. """
    if index.size == 0:
        raise ValueError("Got empty index")
    # Regex pattern "^([ACGT])([0-9]+)$" finds the base and position
    # in each index name.
    pattern = re.compile(f"^([{BASES.decode()}])([0-9]+)$")
    try:
        # Match each index name using the pattern.
        matches = list(map(pattern.match, index))
    except TypeError:
        raise ValueError(f"Invalid indexes: {index}")
    try:
        # Obtain the two groups (base and position) from each match
        # and unzip them into two tuples.
        bases, pos_strs = zip(*map(re.Match.groups, matches))
    except TypeError:
        # TypeError is raised if any match is None, which happens if
        # a column fails to match the pattern.
        invalid_idx = [idx for idx, match in zip(index, matches, strict=True)
                       if match is None]
        raise ValueError(f"Invalid indexes: {invalid_idx}")
    # Join the tuple of bases into a DNA sequence.
    seq = DNA("".join(bases).encode())
    # Cast the tuple of strings of positions into an integer array.
    positions = np.asarray(list(map(int, pos_strs)), dtype=int)
    return seq, positions


def index_to_seq(index: pd.Index):
    return index_to_seq_pos(index)[0]


def index_to_pos(index: pd.Index):
    if index.size == 0:
        # Handle the edge case where index is empty, which would make
        # index_to_seq_pos raise an error.
        return np.array([], dtype=int)
    return index_to_seq_pos(index)[1]


def pos_to_array(pos: Iterable[int]) -> np.ndarray:
    pos = np.asarray(pos)
    # Check for duplicate positions.
    if np.any(dups := np.unique(pos, return_counts=True)[1] > 1):
        raise ValueError(f"Got duplicate positions: {pos[dups]}")
    return pos


class Section(object):
    """ Section of a reference sequence between two coordinates. """

    def __init__(self, /, ref: str, refseq: DNA, *,
                 end5: int | None = None, end3: int | None = None,
                 name: str | None = None):
        """
        Parameters
        ----------
        ref: str
            Name of the reference sequence
        refseq: DNA
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
            Name of the section. If blank, defaults to `self.range`.
        """
        self.ref = ref
        if end5 is None and end3 is None and name is None:
            # Set the name to 'full' if not given.
            name = FULL
        if end5 is None:
            # Set the 5' end to the first position in the reference.
            end5 = 1
        elif end5 < 0:
            # Compute the corresponding positive 5' coordinate.
            end5 += len(refseq) + 1
        if end3 is None:
            # Set the 3' end to the last position in the reference.
            end3 = len(refseq)
        elif end3 < 0:
            # Compute the corresponding positive 3' coordinate.
            end3 += len(refseq) + 1
        if not 1 <= end5 <= end3 <= len(refseq):
            raise ValueError("Must have 1 ≤ end5 ≤ end3 ≤ len(refseq), "
                             f"but got end5 = {end5}, end3 = {end3}, and "
                             f"len(refseq) = {len(refseq)}")
        self.end5: int = end5
        self.end3: int = end3
        self.seq: DNA = refseq[end5 - 1: end3]
        self.full = end5 == 1 and end3 == len(refseq)
        if name is None:
            self.name = self.hyphen
        elif isinstance(name, str):
            self.name = name if name else self.hyphen
        else:
            raise TypeError(f"Parameter 'name' must be 'str', not {type(str)}")
        # Initialize an empty set of masks.
        self._masks: dict[str, np.ndarray] = dict()

    @cached_property
    def range(self):
        """ All numeric positions in the section. """
        return np.arange(self.end5, self.end3 + 1, dtype=int)

    @cached_property
    def range0(self):
        """ All 0-indexed numeric positions in the section. """
        return np.arange(self.length, dtype=int)

    @property
    def focus(self):
        """ Relevant numeric positions in the section. """
        # Alias for unmasked_pos.
        return self.unmasked_pos

    @property
    def focus0(self):
        """ Relevant 0-indexed numeric positions in the section. """
        # Alias for unmasked_pos.
        return self.unmasked_pos - self.end5

    @cached_property
    def range_index(self):
        """ Index names of all positions in the section. """
        return seq_pos_to_index(self.seq, self.range, start=self.end5)

    @property
    def focus_index(self):
        """ Index names of relevant positions in the section. """
        return seq_pos_to_index(self.seq, self.focus, start=self.end5)

    @property
    def length(self):
        """ Length of the entire section. """
        return self.range.size

    @property
    def size(self):
        """ Number of relevant positions in the section. """
        return self.length - self.masked_pos.size

    @property
    def coord(self):
        """ Tuple of the 5' and 3' coordinates. """
        return self.end5, self.end3

    @property
    def hyphen(self):
        return f"{self.end5}-{self.end3}"

    @property
    def ref_sect(self):
        return f"{self.ref}:{self.name}"

    def to_dict(self):
        return {"ref": self.ref, "seq": self.seq, "sect": self.name,
                "end5": self.end5, "end3": self.end3}

    MASK_POLYA = "pos-polya"
    MASK_GU = "pos-gu"
    MASK_POS = "pos-user"

    @property
    def mask_names(self):
        """ Names of the masks. """
        return list(self._masks)

    @property
    def masked_pos(self) -> np.ndarray:
        """ Masked positions as integers. """
        # Do not cache this method since self._masks can change.
        return reduce(np.union1d, self._masks.values(), np.array([], dtype=int))

    @property
    def masked_bool(self) -> np.ndarray:
        """ Masked positions as a boolean array. """
        # Do not cache this method since self.masked_int can change.
        return np.isin(self.range, self.masked_pos)

    @property
    def unmasked_bool(self) -> np.ndarray:
        """ Unmasked positions as a boolean array. """
        # Do not cache this method since self.masked_bool can change.
        return np.logical_not(self.masked_bool)

    @property
    def unmasked_pos(self) -> np.ndarray:
        """ Unmasked positions as integers. """
        # Do not cache this method since self.unmasked_bool can change.
        return self.range[self.unmasked_bool]

    def get_mask(self, name: str):
        """ Get the positions masked under the given name. """
        return self._masks[name]

    def add_mask(self, name: str, mask_pos: Iterable[int]):
        """ Mask the integer positions in the array `mask_pos`. """
        if name in self._masks:
            raise ValueError(f"Mask '{name}' was already set")
        # Convert positions to a NumPy integer array.
        pos = np.asarray(mask_pos, dtype=int)
        # Check for positions outside the section.
        if np.any(pos < self.end5) or np.any(pos > self.end3):
            invalid = pos[np.logical_or(pos < self.end5, pos > self.end3)]
            logger.warning(f"Got positions to mask ouside of {self}: {invalid}")
        # Record the positions that have not already been masked.
        self._masks[name] = pos[np.isin(pos, self.masked_pos, invert=True)]

    def _find_gu(self, exclude_gu: bool = True) -> np.ndarray:
        """ Array of each position whose base is neither A nor C. """
        if exclude_gu:
            # Convert the sequence to an array for broadcasting.
            seq_array = self.seq.to_int_array()
            # Mark whether each position is neither A nor C.
            gu_pos = np.logical_and(seq_array != A_INT, seq_array != C_INT)
            # Return the integer positions.
            return self.range[gu_pos]
        # Mask no positions.
        return np.array([], dtype=int)

    def mask_gu(self, exclude: bool):
        """ Mask positions whose base is neither A nor C. """
        self.add_mask(self.MASK_GU, self._find_gu(exclude))

    def _find_polya(self, min_length: int) -> np.ndarray:
        """ Array of each position within a stretch of `min_length` or
        more consecutive adenines. """
        if min_length < 0:
            raise ValueError(f"min_length must be ≥ 0, but got {min_length}")
        # Initialize a list of 0-indexed positions in poly(A) sequences.
        polya_pos = list()
        if min_length > 0:
            # Generate a pattern that matches stretches of consecutive
            # adenines that are at least as long as min_length.
            polya_pattern = b"%c{%i,}" % (A_INT, min_length)
            # Add the 0-indexed positions in every poly(A) sequence.
            for polya in re.finditer(polya_pattern, self.seq):
                polya_pos.extend(range(polya.start(), polya.end()))
        # Convert the positions to an array with natural indexing.
        return np.array(polya_pos, dtype=int) + self.end5

    def mask_polya(self, min_length: int):
        """ Mask poly(A) stretches with length ≥ `min_length`. """
        self.add_mask(self.MASK_POLYA, self._find_polya(min_length))

    def mask_pos(self, pos: Iterable[int]):
        """ Mask arbitrary positions. """
        self.add_mask(self.MASK_POS, pos)

    def __str__(self):
        return f"Section {self.ref_sect} ({self.hyphen}) {self.mask_names}"

    def __eq__(self, other):
        if self is other:
            # If self and other are the same object, they must be equal.
            return True
        if not isinstance(other, Section):
            # Cannot compare to an object that is not a Section.
            return NotImplemented
        # Compare the sections' sequences, positions, and names.
        if not (self.ref == other.ref and self.seq == other.seq
                and self.end5 == other.end5 and self.end3 == other.end3
                and self.name == other.name):
            return False
        # If that comparison passed, then compare their mask names.
        if sorted(self.mask_names) != sorted(other.mask_names):
            return False
        # Compare all mask values and return False if any differ.
        for name, mask in self._masks.items():
            if not np.array_equal(mask, other.get_mask(name)):
                return False
        # All checks for equality passed.
        return True


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
       - (primer_gap + 1) if rev is given, else the length of refseq
    """

    def __init__(self, /, ref: str, refseq: DNA | None, *,
                 name: str | None = None,
                 end5: int | None = None, end3: int | None = None,
                 fwd: DNA | None = None, rev: DNA | None = None,
                 primer_gap: int | None = None):
        """
        Parameters
        ----------
        ref: str
            see superclass
        refseq: DNA
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
        if refseq is None:
            raise ValueError(f"No sequence for reference named '{ref}'. Check "
                             "that you gave the right reference sequence file "
                             "and spelled the name of the reference correctly.")
        if end5 is None:
            # No 5' end coordinate was given.
            if fwd is not None:
                # A forward primer was given.
                if primer_gap is None:
                    raise TypeError("Must give primer_gap if using primers")
                # Locate the end of the forward primer, then start the
                # section (primer_gap + 1) positions downstream.
                end5 = self.locate(refseq, fwd).pos3 + (primer_gap + 1)
        if end3 is None:
            # No 3' end coordinate was given.
            if rev is not None:
                # A reverse primer was given.
                if primer_gap is None:
                    raise TypeError("Must give primer_gap if using primers")
                # Locate the start of the reverse primer, then end the
                # section (primer_gap + 1) positions upstream.
                end3 = self.locate(refseq, rev.rc).pos5 - (primer_gap + 1)
        super().__init__(refseq=refseq, ref=ref,
                         end5=end5, end3=end3, name=name)

    @staticmethod
    def locate(refseq: DNA, primer: DNA) -> SectionTuple:
        """
        Return the 5' and 3' positions (1-indexed) of a primer within a
        reference sequence. The primer must occur exactly once in the
        reference, otherwise an error is raised.

        Parameters
        ----------
        refseq: DNA
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
        matches = list(re.finditer(primer, refseq))
        if not matches:
            raise ValueError(f"Primer '{primer}' is not in ref '{refseq}'")
        if len(matches) > 1:
            raise ValueError(f"Primer '{primer}' occurs {len(matches)} times "
                             f"in ref '{refseq}'")
        # Add 1 to convert from 0-indexed to 1-indexed.
        pos5 = matches[0].start() + 1
        # No change is needed to convert from exclusive 0-indexed to
        # inclusive 1-indexed.
        pos3 = matches[0].end()
        return SectionTuple(pos5, pos3)


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
                 refseqs: Iterable[tuple[str, DNA]], *,
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
        for ref, seq in refseqs:
            self._sections[ref] = dict()
            for end5, end3 in ref_coords[ref]:
                # Add a section for each pair of 5' and 3' coordinates.
                self._add_section(ref=ref, refseq=seq, end5=end5, end3=end3,
                                  primer_gap=primer_gap,
                                  name=section_names.get((ref, end5, end3)))
            for fwd, rev in ref_primers[ref]:
                # Add a section for each pair of fwd and rev primers.
                self._add_section(ref=ref, refseq=seq, fwd=fwd, rev=rev,
                                  primer_gap=primer_gap)
            if not self._sections[ref]:
                # If no sections were given for the reference, then add
                # a section that spans the full reference.
                self._add_section(ref=ref, refseq=seq,
                                  primer_gap=primer_gap)
        if extra := (set(ref_coords) | set(ref_primers)) - set(self._sections):
            logger.warning(f"No sequences given for references {sorted(extra)}")

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
