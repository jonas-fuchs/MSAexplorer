"""
this contains the dataclasses used to store the data for the msa explorer. these are not meant to be used outside of this package.
"""

# build-in
from dataclasses import dataclass

# libs
from numpy import ndarray

@dataclass(frozen=True)
class PairwiseDistanceResult:
    """
    Result container for pairwise identity values between a reference/consensus
    sequence and each sequence in the alignment.

    The object remains iterable so it can be unpacked as a tuple:
    ``reference_label, sequence_ids, distances = result``.
    """

    reference_id: str
    sequence_ids: list[str]
    distances: ndarray

    def __iter__(self):
        yield self.reference_id
        yield self.sequence_ids
        yield self.distances