"""
Shared pytest fixtures and helper functions for MSAexplorer tests.
"""

import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment


def create_alignment(sequences_dict):
    """
    Helper function to create a MultipleSeqAlignment from a dictionary.

    Args:
        sequences_dict: Dictionary with sequence IDs as keys and sequences as values

    Returns:
        Bio.Align.MultipleSeqAlignment object
    """
    records = [SeqRecord(Seq(seq), id=seq_id) for seq_id, seq in sequences_dict.items()]
    return MultipleSeqAlignment(records)

