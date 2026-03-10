"""
Tests for consensus sequence generation.
"""

import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from msaexplorer.explore import MSA


def create_alignment(sequences_dict):
    """Helper function to create a MultipleSeqAlignment from a dictionary."""
    records = [SeqRecord(Seq(seq), id=seq_id) for seq_id, seq in sequences_dict.items()]
    return MultipleSeqAlignment(records)


class TestGetConsensus:
    """Test the get_consensus method with various parameters."""

    def test_consensus_no_threshold_dna(self):
        """Test consensus with no threshold returns most frequent nucleotide at each position."""
        alignment_dict = {
            'seq1': 'ATGC',
            'seq2': 'ATGC',
            'seq3': 'ATAC',
            'seq4': 'AAGC'
        }
        aln = MSA(create_alignment(alignment_dict))
        consensus = aln.get_consensus()

        assert consensus == 'ATGC'

    def test_consensus_with_threshold_no_ambig(self):
        """Test consensus with threshold but without ambiguous characters."""
        alignment_dict = {
            'seq1': 'ATGC',
            'seq2': 'AAGC',
            'seq3': 'ACGC',
            'seq4': 'ATGC'
        }
        aln = MSA(create_alignment(alignment_dict))
        consensus = aln.get_consensus(threshold=0.6)

        assert consensus[1] == 'N'  # No single nucleotide reaches 60%


    def test_consensus_with_threshold_and_ambig_dna(self):
        """Test consensus with threshold and ambiguous characters for DNA."""
        alignment_dict = {
            'seq1': 'ATGC',
            'seq2': 'AAGC',
            'seq3': 'ACGC',
            'seq4': 'ATGC'
        }
        aln = MSA(create_alignment(alignment_dict))
        consensus = aln.get_consensus(threshold=0.8, use_ambig_nt=True)

        assert consensus[1] == 'H'  # A, C, T


    def test_consensus_simple_ambig_dna(self):
        """Test consensus with simple two-way ambiguity in DNA."""
        alignment_dict = {
            'seq1': 'AGGG',
            'seq2': 'GGGG',
            'seq3': 'AGGG',
            'seq4': 'GGGG'
        }
        aln = MSA(create_alignment(alignment_dict))
        consensus = aln.get_consensus(threshold=0.9, use_ambig_nt=True)

        assert consensus[0] == 'R'  # A and G -> R

    def test_consensus_threshold_zero(self):
        """Test consensus with threshold of 0 returns most frequent base."""
        alignment_dict = {
            'seq1': 'ATGC',
            'seq2': 'AAGC',
            'seq3': 'GGGG'
        }
        aln = MSA(create_alignment(alignment_dict))
        consensus = aln.get_consensus(threshold=0)

        assert consensus[1] in ['A', 'T', 'G']  # All equal frequency


    def test_consensus_with_gaps(self):
        """Test consensus handles gaps correctly."""
        alignment_dict = {
            'seq1': 'A-GC',
            'seq2': 'ATGC',
            'seq3': 'ATGC',
            'seq4': 'AT-C'
        }
        aln = MSA(create_alignment(alignment_dict))
        consensus = aln.get_consensus(threshold=0.5)

        assert consensus[1] == 'T'  # T=75% reaches 50%
        assert consensus[2] == 'G'  # G=75% reaches 50%

    def test_consensus_with_ambig_in_sequences(self):
        """Test consensus handles ambiguous characters in input sequences."""
        alignment_dict = {
            'seq1': 'ATGC',
            'seq2': 'ARGC',  # R = A or G
            'seq3': 'ATGC',
            'seq4': 'ATGC'
        }
        aln = MSA(create_alignment(alignment_dict))
        consensus = aln.get_consensus()

        assert consensus[1] == 'T'

    def test_consensus_rna(self):
        """Test consensus with RNA sequences."""
        alignment_dict = {
            'seq1': 'AUGC',
            'seq2': 'AUGC',
            'seq3': 'AUAC',
            'seq4': 'AAGC'
        }
        aln = MSA(create_alignment(alignment_dict))
        consensus = aln.get_consensus()

        assert consensus == 'AUGC'

    def test_consensus_rna_with_ambig(self):
        """Test consensus with RNA and ambiguous characters."""
        alignment_dict = {
            'seq1': 'AGGG',
            'seq2': 'UGGG',
            'seq3': 'AGGG',
            'seq4': 'UGGG'
        }
        aln = MSA(create_alignment(alignment_dict))
        consensus = aln.get_consensus(threshold=0.9, use_ambig_nt=True)

        assert consensus[0] == 'W'

    def test_consensus_amino_acid(self):
        """Test consensus with amino acid sequences."""
        alignment_dict = {
            'seq1': 'ARND',
            'seq2': 'ARND',
            'seq3': 'ACND',
            'seq4': 'ARND'
        }
        aln = MSA(create_alignment(alignment_dict))
        consensus = aln.get_consensus(threshold=0.6)

        assert consensus[1] == 'R'  # R=75% reaches 60%


    def test_consensus_amino_acid_with_x(self):
        """Test consensus with amino acid sequences using X for ambiguity."""
        alignment_dict = {
            'seq1': 'ARND',
            'seq2': 'ACND',
            'seq3': 'AGND',
            'seq4': 'AKND'
        }
        aln = MSA(create_alignment(alignment_dict))
        # Position 1: all different amino acids, none reach 60%
        consensus = aln.get_consensus(threshold=0.6)

        assert consensus[1] == 'X'  # No single AA reaches 60%

    def test_consensus_threshold_validation(self):
        """Test that invalid threshold values raise ValueError."""
        alignment_dict = {
            'seq1': 'ATGC',
            'seq2': 'ATGC'
        }
        aln = MSA(create_alignment(alignment_dict))

        with pytest.raises(ValueError, match='Threshold must be between 0 and 1'):
            aln.get_consensus(threshold=-0.1)

        with pytest.raises(ValueError, match='Threshold must be between 0 and 1'):
            aln.get_consensus(threshold=1.5)

    def test_consensus_ambig_without_threshold_raises_error(self):
        """Test that using ambig characters without threshold raises ValueError."""
        alignment_dict = {
            'seq1': 'ATGC',
            'seq2': 'ATGC'
        }
        aln = MSA(create_alignment(alignment_dict))

        with pytest.raises(ValueError, match='To calculate ambiguous nucleotides, set a threshold > 0'):
            aln.get_consensus(use_ambig_nt=True)

    def test_consensus_ambig_with_amino_acid_raises_error(self):
        """Test that using ambig characters with amino acids raises ValueError."""
        alignment_dict = {
            'seq1': 'ARND',
            'seq2': 'ARND'
        }
        aln = MSA(create_alignment(alignment_dict))

        with pytest.raises(ValueError, match='Ambiguous characters can not be calculated for amino acid alignments'):
            aln.get_consensus(threshold=0.5, use_ambig_nt=True)

