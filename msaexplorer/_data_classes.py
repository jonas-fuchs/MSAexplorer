"""
this contains the dataclasses used to store the data for the msa explorer. these are not meant to be used outside of this package.
"""

# build-in
from dataclasses import dataclass

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
    aln_type: str
    reference_id: str | None

    def __post_init__(self):
        if self.positions.shape != self.values.shape:
            raise ValueError("positions and values must have the same shape")

    # dunder methods
    def __len__(self) -> int:
        return len(self.values)

    def __getitem__(self, index: int) -> float:
        return self.values[index]

    def __iter__(self):
        yield self.stat_name
        yield self.positions
        yield self.values

    def __contains__(self, item: float) -> bool:
        return item in self.positions

    # normal mehods
    def as_array(self) -> ndarray:
        return self.values

    def as_list(self) -> list:
        return self.values.tolist()


@dataclass(frozen=True)
class PairwiseDistance:
    """
    Result container for Pairwise distances. Array can either be a 2D (compared to reference) or 3D array.
    """

    reference_id: str | None
    sequence_ids: list[str]
    distances: ndarray

    # dunder methods
    def __iter__(self):
        yield self.reference_id
        yield self.sequence_ids
        yield self.distances

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
