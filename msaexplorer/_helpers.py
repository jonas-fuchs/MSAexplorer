"""
This contains helper functions, not intended to be used outside of this package.
"""

import os, io, math
import numpy as np
from msaexplorer import config
from typing import Callable, Dict
from matplotlib.colors import is_color_like
from Bio.Align import MultipleSeqAlignment
from Bio import AlignIO


# explore helpers
def _get_line_iterator(source):
    """
    allow reading in both raw string or paths
    """
    if isinstance(source, str) and os.path.exists(source):
        return open(source, 'r')
    else:
        return io.StringIO(source)


def _read_alignment(source: str | MultipleSeqAlignment) -> dict:
    """
    Parse MSA alignment using Biopython with automatic format detection.
    :param source: Path to alignment file, raw alignment string, or Bio.Align.MultipleSeqAlignment object
    :possible_chars: list of possible characters in the alignment
    :return: dictionary with ids as keys and sequences as values
    """
    # Handle Bio.Align.MultipleSeqAlignment objects
    if isinstance(source, MultipleSeqAlignment):
        aln_dict = {record.id: str(record.seq).upper() for record in source}
    else:
        # Try multiple formats in order of likelihood
        formats_to_try = ["fasta", "clustal", "phylip", "stockholm", "nexus"]

        aln_dict = None
        for fmt in formats_to_try:
            try:
                with _get_line_iterator(source) as handle:
                    alignment = AlignIO.read(handle, fmt)
                    aln_dict = {record.id: str(record.seq).upper() for record in alignment}
                    break
            except Exception:
                continue

        if aln_dict is None:
            # If no format worked, raise an error
            raise ValueError(f"Alignment file could not be parsed. Supported formats: {', '.join(formats_to_try)}")

    # Validate alignment
    if not aln_dict:
        raise ValueError(f"Alignment does not contain any sequences.")

    if len(aln_dict) < 2:
        raise ValueError("Alignment must contain more than one sequence.")

    # Check for non-allowed characters
    for sequence_id, seq in aln_dict.items():
        invalid_chars = set(seq) - set(config.POSSIBLE_CHARS)
        if invalid_chars:
            raise ValueError(
                f"{sequence_id} contains invalid characters: {', '.join(invalid_chars)}. Allowed chars are: {config.POSSIBLE_CHARS}."
            )

    # Validate all sequences have same length
    first_seq_len = len(next(iter(aln_dict.values())))
    for sequence_id, seq in aln_dict.items():
        if len(seq) != first_seq_len:
            raise ValueError(
                f"All alignment sequences must have the same length. Sequence '{sequence_id}' has length {len(seq)}, expected {first_seq_len}."
            )

    return aln_dict


def _create_distance_calculation_function_mapping() -> Dict[str, Callable[[str, str, int], float]]:
    """
    create a mapping of distance types to distance calculation functions
    :return: dictionary of mappings
    """

    def ghd(seq1: str, seq2: str, aln_length: int) -> float:
        """
        global hamming distance - defined as percentage of total number of matches
        """
        return sum(c1 == c2 for c1, c2 in zip(seq1, seq2)) / aln_length * 100

    def lhd(seq1: str, seq2: str, aln_length: int) -> float:
        """
        local hamming distance - defined as total number of matches excluding terminal gaps
        """
        # Trim gaps from both sides
        i, j = 0, aln_length - 1
        while i < aln_length and (seq1[i] == '-' or seq2[i] == '-'):
            i += 1
        while j >= 0 and (seq1[j] == '-' or seq2[j] == '-'):
            j -= 1
        if i > j:
            return 0.0

        seq1_, seq2_ = seq1[i:j + 1], seq2[i:j + 1]
        matches = sum(c1 == c2 for c1, c2 in zip(seq1_, seq2_))
        length = j - i + 1

        return (matches / length) * 100 if length > 0 else np.nan

    def ged(seq1: str, seq2: str, aln_length: int = None) -> float:
        """
        gap excluded distance - defined as percentage of total number of matches excluding all gaps
        """

        diff, total = 0, 0

        for c1, c2 in zip(seq1, seq2):
            if c1 != '-' and c2 != '-':
                total += 1
                if c1 != c2:
                    diff += 1

        return (1 - diff / total) * 100 if total > 0 else np.nan

    def gcd(seq1: str, seq2: str, aln_length: int = None) -> float:
        """
        gap compressed distance - defined as percentage of total number of matches with sequential gap mismatches
        counting as a single mismatch
        """
        matches = 0
        mismatches = 0
        in_gap = False

        for char1, char2 in zip(seq1, seq2):
            if char1 == '-' and char2 == '-':  # Shared gap: do nothing
                continue
            elif char1 == '-' or char2 == '-':  # Gap in only one sequence
                if not in_gap:  # Start of a new gap stretch
                    mismatches += 1
                    in_gap = True
            else:  # No gaps
                in_gap = False
                if char1 == char2:  # Matching characters
                    matches += 1
                else:  # Mismatched characters
                    mismatches += 1

        return matches / (matches + mismatches) * 100 if (matches + mismatches) > 0 else np.nan

    def jc69(seq1: str, seq2: str, aln_length: int = None) -> float:
        """
        Jukes-Cantor 1969 (JC69) corrected identity.
        Gaps are excluded. The proportion of differing sites (p-distance) is corrected
        for multiple hits: d = -(3/4) * ln(1 - (4/3) * p).
        Returns (1 - d) * 100 as a corrected percent identity (100 = identical).
        Returns 0 when p >= 0.75 (formula undefined / sequence saturated).
        """
        diff, total = 0, 0
        for c1, c2 in zip(seq1, seq2):
            if c1 != '-' and c2 != '-':
                total += 1
                if c1 != c2:
                    diff += 1
        if total == 0:
            return np.nan
        p = diff / total
        if p == 0.0:
            return 100.0
        correction = 1.0 - (4.0 / 3.0) * p
        if correction <= 0.0:  # saturated – formula undefined
            return 0.0
        d = -(3.0 / 4.0) * math.log(correction)
        return max(0.0, (1.0 - d) * 100.0)

    def k2p(seq1: str, seq2: str, aln_length: int = None) -> float:
        """
        Kimura 2-Parameter (K2P / K80) corrected identity.
        Gaps are excluded. Transitions and transversions (Tv) are weighted separately:
        d = -(1/2) * ln(1 - 2P - Q) - (1/4) * ln(1 - 2Q)
        where P = Ti / total and Q = Tv / total.
        Returns (1 - d) * 100 as a corrected percent identity (100 = identical).
        Returns 0 when the logarithm arguments become non-positive (saturated).
        """
        transitions = [{'A', 'G'}, {'C', 'T'}]

        ts, tv, total = 0, 0, 0
        for c1, c2 in zip(seq1, seq2):
            if c1 != '-' and c2 != '-':
                total += 1
                if c1 != c2:
                    if {c1, c2} in transitions:
                        ts += 1
                    else:
                        tv += 1
        if total == 0:
            return np.nan
        if ts == 0 and tv == 0:
            return 100.0
        P = ts / total   # transition proportion
        Q = tv / total   # transversion proportion
        term1 = 1.0 - 2.0 * P - Q
        term2 = 1.0 - 2.0 * Q
        # saturated – formula undefined
        if term1 <= 0.0 or term2 <= 0.0:
            return 0.0
        # calculate distance
        d = -0.5 * math.log(term1) - 0.25 * math.log(term2)

        return max(0.0, (1 - d) * 100.0)

    # Map distance type to corresponding function
    distance_functions: Dict[str, Callable[[str, str, int], float]] = {
        'ghd': ghd,
        'lhd': lhd,
        'ged': ged,
        'gcd': gcd,
        'jc69': jc69,
        'k2p': k2p,
    }

    return distance_functions

# export helpers
def _check_and_create_path(path: str):
    """
    Check and create path if it doesn't exist.
    :param path: string to file
    """
    if path is not None:
        output_dir = os.path.dirname(path)
        if output_dir and not os.path.exists(output_dir):
            os.makedirs(output_dir)

# draw helpers
def _validate_color(c):
    """
    validate color and raise error
    """
    if not is_color_like(c):
        raise ValueError(f'{c} is not a color')

def _get_contrast_text_color(rgba_color):
    """
    compute the brightness of a color
    """
    r, g, b, a = rgba_color
    brightness = (r * 299 + g * 587 + b * 114) / 1000

    return 'white' if brightness < 0.5 else 'black'

