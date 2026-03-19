"""
this contains the dataclasses used to store the data for the msa explorer. these are not meant to be used outside of this package.
"""

# build-in
from dataclasses import dataclass, field

# libs
from numpy import ndarray


@dataclass(frozen=True)
class AlignmentStats:
    """
    Generic result container for position-based statistics.
    """

    stat_name: str
    positions: ndarray
    values: ndarray

    def __post_init__(self):
        if self.positions.shape != self.values.shape:
            raise ValueError("positions and values must have the same shape")

    # dunder methods
    def __len__(self) -> int:
        return len(self.values)

    def __getitem__(self, index: int) -> float:
        return self.values[index]

    def __contains__(self, item: float) -> bool:
        return item in self.positions


@dataclass(frozen=True)
class PairwiseDistance:
    """
    Result container for Pairwise distances. Array can either be a 2D (compared to reference) or 3D array.
    """

    reference_id: str | None
    sequence_ids: list[str]
    distances: ndarray

    # dunder methods
    def __len__(self) -> int:
        return len(self.sequence_ids)

    def __getitem__(self, index: int | str) -> float:
        """
        Different ways to access the distance matrix.
        - pd[0] -> Distance first sequence to all other sequences
        - pd['seq_name'] -> Distance of a specific sequence to all other sequences
        """
        if isinstance(index, str):
            idx = self.sequence_ids.index(index)
            return self.distances[idx]
        return self.distances[index]

    def __contains__(self, item: str) -> bool:
        """Seq ID present"""
        return item in self.sequence_ids


@dataclass(frozen=True)
class OpenReadingFrame:
    """
    Represents a single conserved ORF detected across an alignment.

    Attributes:
        orf_id:       Unique identifier, e.g. ``'ORF_0'``.
        location:     Main ORF boundaries as a tuple of ``(start, stop)`` pairs
                      (0-based, half-open).  Typically a single pair, but may
                      carry additional coordinates for split ORFs.
        frame:        Reading frame (0, 1, or 2).
        strand:       ``'+'`` for forward, ``'-'`` for reverse complement.
        conservation: Percentage of fully identical alignment columns inside the ORF.
        internal:     Tuple of ``(start, stop)`` pairs for nested (internal) ORFs
                      that share the same stop codon.
    """

    orf_id: str
    location: tuple[tuple[int, int], ...]
    frame: int
    strand: str
    conservation: float
    internal: tuple[tuple[int, int], ...] = field(default_factory=tuple)

    def __post_init__(self):
        if self.strand not in ('+', '-'):
            raise ValueError(f"strand must be '+' or '-', got {self.strand!r}")
        if not (0 <= self.frame <= 2):
            raise ValueError(f"frame must be 0, 1, or 2, got {self.frame!r}")

    def __len__(self) -> int:
        """Length of the main ORF in alignment columns."""
        return self.location[0][1] - self.location[0][0]

    def __contains__(self, position: int) -> bool:
        """True if *position* (0-based) falls inside the main ORF."""
        start, stop = self.location[0]
        return start <= position < stop


@dataclass(frozen=True)
class OrfContainer:
    """
    Ordered collection of `OpenReadingFrame` objects returned by
    `MSA.get_conserved_orfs` or`MSA.get_non_overlapping_conserved_orfs`.

    The class intentionally mimics a *read-only dict* interface
    """

    orfs: tuple[OpenReadingFrame, ...] = field(default_factory=tuple)

    # dict-like interface
    def keys(self) -> list[str]:
        """Return ORF identifiers in insertion order."""
        return [orf.orf_id for orf in self.orfs]

    def values(self) -> list[OpenReadingFrame]:
        """Return :class:`OpenReadingFrame` objects in insertion order."""
        return list(self.orfs)

    def items(self) -> list[tuple[str, OpenReadingFrame]]:
        """Return ``(orf_id, OpenReadingFrame)`` pairs in insertion order."""
        return [(orf.orf_id, orf) for orf in self.orfs]

    # dunder methods
    def __len__(self) -> int:
        return len(self.orfs)

    def __bool__(self) -> bool:
        return len(self.orfs) > 0

    def __iter__(self):
        """Iterate over ORF identifiers (mirrors ``dict.__iter__``)."""
        return iter(orf.orf_id for orf in self.orfs)

    def __getitem__(self, key: str | int) -> OpenReadingFrame:
        """
        Access an ORF by identifier string or integer index.

        Examples::
            orfs['ORF_0']   # by identifier
            orfs[0]         # by index
        """
        if isinstance(key, int):
            return self.orfs[key]
        for orf in self.orfs:
            if orf.orf_id == key:
                return orf
        raise KeyError(key)

    def __contains__(self, orf_id: str) -> bool:
        """True if an ORF with *orf_id* is present in the collection."""
        return any(orf.orf_id == orf_id for orf in self.orfs)

