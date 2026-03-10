"""
Tests for conserved ORF detection.
"""

import pytest
from msaexplorer.explore import MSA
from conftest import create_alignment


class TestGetConservedOrfs:
    """Test the get_conserved_orfs method."""

    def test_basic_orf_detection(self):
        """Test detection of a simple conserved ORF."""
        # ATG = start, TAA = stop
        alignment_dict = {
            'seq1': 'ATGAAAAAATAA',
            'seq2': 'ATGAAAAAATAA',
            'seq3': 'ATGAAAAAATAA'
        }
        aln = MSA(create_alignment(alignment_dict))
        orfs = aln.get_conserved_orfs(min_length=9)

        assert len(orfs) > 0
        orf_key = list(orfs.keys())[0]
        assert orfs[orf_key]['strand'] == '+'
        assert orfs[orf_key]['location'][0][0] == 0  # Start at position 0
        assert orfs[orf_key]['location'][0][1] == 12  # End at position 9

    def test_orf_with_gaps(self):
        """Test ORF detection handles gaps correctly."""
        alignment_dict = {
            'seq1': 'ATG---GGGAAATAA',
            'seq2': 'ATGAAAGGG---TAA',
            'seq3': 'ATGAAAGGG---TAA'
        }
        aln = MSA(create_alignment(alignment_dict))
        orfs = aln.get_conserved_orfs(min_length=9)

        assert len(orfs) == 1

    def test_no_orfs(self):
        """Test that empty dict is returned when no ORFs found."""
        alignment_dict = {
            'seq1': 'AAAAAAAAAAA',
            'seq2': 'AAAAAAAAAAA',
            'seq3': 'AAAAAAAAAAA'
        }
        aln = MSA(create_alignment(alignment_dict))
        orfs = aln.get_conserved_orfs(min_length=9)

        assert len(orfs) == 0

    def test_orf_with_internal_stop(self):
        """Test that ORFs with internal stop codons are rejected."""
        alignment_dict = {
            'seq1': 'ATGTAAATGTAA',
            'seq2': 'ATGTAAATGTAA',
            'seq3': 'ATGTAAATGTAA'
        }
        aln = MSA(create_alignment(alignment_dict))
        orfs = aln.get_conserved_orfs(min_length=9)

        assert len(orfs) == 0  # Should not find ORF because of internal stop at position 3-5

    def test_orf_multiple_frames(self):
        """Test ORF detection across different reading frames."""
        # One ORF in frame 1
        alignment_dict = {
            'seq1': 'CATGAAATAAGATGCCCTAG',
            'seq2': 'CATGAAATAAGATGCCCTAG',
            'seq3': 'CATGAAATAAGATGCCCTAG'
        }
        aln = MSA(create_alignment(alignment_dict))
        orfs = aln.get_conserved_orfs(min_length=9)

        # Should detect ORFs in frame 1
        frames = [orf['frame'] for orf in orfs.values()]
        assert 1 in frames

    def test_orf_reverse_strand(self):
        """Test ORF detection on reverse complement strand."""
        # Should find ORF in both strands
        alignment_dict = {
            'seq1': 'ATGTTATTTCATTAA',
            'seq2': 'ATGTTATTTCATTAA',
            'seq3': 'ATGTTATTTCATTAA'
        }
        aln = MSA(create_alignment(alignment_dict))
        orfs = aln.get_conserved_orfs(min_length=9)
        strands = [orf['strand'] for orf in orfs.values()]

        assert strands == ['+', '-']

    def test_orf_non_conserved_start(self):
        """Test that non-conserved start codons are rejected."""
        alignment_dict = {
            'seq1': 'ATGAAATAA',  # Has ATG
            'seq2': 'ATGAAATAA',  # Has ATG
            'seq3': 'TTGAAATAA'   # Has TTG (not a start codon)
        }
        aln = MSA(create_alignment(alignment_dict))
        orfs = aln.get_conserved_orfs(min_length=9)

        assert len(orfs) == 0

    def test_orf_non_conserved_stop(self):
        """Test that non-conserved stop codons are rejected."""
        alignment_dict = {
            'seq1': 'ATGAAATAA',  # Has TAA stop
            'seq2': 'ATGAAATAA',  # Has TAA stop
            'seq3': 'ATGAAACAA'   # Has CAA (not a stop)
        }
        aln = MSA(create_alignment(alignment_dict))
        orfs = aln.get_conserved_orfs(min_length=9)

        assert len(orfs) == 0

    def test_orf_identity_cutoff(self):
        """Test ORF filtering by identity/conservation cutoff."""
        alignment_dict = {
            'seq1': 'ATGAAATAA',
            'seq2': 'ATGAAATAA',
            'seq3': 'ATGCCCTAA'  # Different middle region
        }
        aln = MSA(create_alignment(alignment_dict))

        # With low cutoff
        orfs_low = aln.get_conserved_orfs(min_length=9, identity_cutoff=50.0)
        assert len(orfs_low) == 1

        # With high cutoff
        orfs_high = aln.get_conserved_orfs(min_length=9, identity_cutoff=99.0)
        assert len(orfs_high) == 0

    def test_orf_internal_orfs_detected(self):
        """Test that internal ORFs are properly marked."""
        # Two start codons to same stop --> one internal
        alignment_dict = {
            'seq1': 'ATGAAAATGCCCTAA',
            'seq2': 'ATGAAAATGCCCTAA',
            'seq3': 'ATGAAAATGCCCTAA'
        }
        aln = MSA(create_alignment(alignment_dict))
        orfs = aln.get_conserved_orfs(min_length=9)

        assert len(orfs['ORF_0']['internal']) == 1


    def test_orf_rna_alignment(self):
        """Test ORF detection works with RNA alignments (U instead of T)."""
        alignment_dict = {
            'seq1': 'AUGAAAUAA',
            'seq2': 'AUGAAAUAA',
            'seq3': 'AUGAAAUAA'
        }
        aln = MSA(create_alignment(alignment_dict))
        orfs = aln.get_conserved_orfs(min_length=9)

        assert len(orfs) == 1

    def test_orf_amino_acid_raises_error(self):
        """Test that amino acid alignments raise TypeError."""
        alignment_dict = {
            'seq1': 'ARNDCQEGHILKMFPSTWYV',
            'seq2': 'ARNDCQEGHILKMFPSTWYV',
            'seq3': 'ARNDCQEGHILKMFPSTWYV'
        }
        aln = MSA(create_alignment(alignment_dict))

        with pytest.raises(TypeError, match='ORF search only for RNA/DNA alignments'):
            aln.get_conserved_orfs(min_length=9)

    def test_orf_invalid_min_length(self):
        """Test that invalid min_length values raise ValueError."""
        alignment_dict = {
            'seq1': 'ATGAAATAA',
            'seq2': 'ATGAAATAA',
            'seq3': 'ATGAAATAA'
        }
        aln = MSA(create_alignment(alignment_dict))

        # Too small
        with pytest.raises(ValueError, match='min_length must be between 6 and'):
            aln.get_conserved_orfs(min_length=3)

        # Larger than the alignment length
        with pytest.raises(ValueError, match='min_length must be between 6 and'):
            aln.get_conserved_orfs(min_length=1000)

    def test_orf_invalid_identity_cutoff(self):
        """Test that invalid identity_cutoff values raise ValueError."""
        alignment_dict = {
            'seq1': 'ATGAAATAA',
            'seq2': 'ATGAAATAA',
            'seq3': 'ATGAAATAA'
        }
        aln = MSA(create_alignment(alignment_dict))

        # Negative cutoff
        with pytest.raises(ValueError, match='conservation cutoff must be between 0 and 100'):
            aln.get_conserved_orfs(min_length=9, identity_cutoff=-1.0)

        # > 100 cutoff
        with pytest.raises(ValueError, match='conservation cutoff must be between 0 and 100'):
            aln.get_conserved_orfs(min_length=9, identity_cutoff=101.0)

    def test_orf_return_structure(self):
        """Test that returned ORF dictionary has correct structure."""
        alignment_dict = {
            'seq1': 'ATGAAATAA',
            'seq2': 'ATGAAATAA',
            'seq3': 'ATGAAATAA'
        }
        aln = MSA(create_alignment(alignment_dict))
        orfs = aln.get_conserved_orfs(min_length=9)

        assert isinstance(orfs, dict)

        # first orf
        orf = orfs[list(orfs.keys())[0]]

        # Check required keys
        assert 'location' in orf
        assert 'frame' in orf
        assert 'strand' in orf
        assert 'conservation' in orf
        assert 'internal' in orf

        # Check types
        assert isinstance(orf['location'], list)
        assert isinstance(orf['frame'], int)
        assert orf['strand'] in ['+', '-']
        assert isinstance(orf['conservation'], float)
        assert isinstance(orf['internal'], list)


class TestGetNonOverlappingConservedOrfs:
    """Test the get_non_overlapping_conserved_orfs method."""

    def test_non_overlapping_basic(self):
        """Test basic non-overlapping ORF selection."""
        # Two non-overlapping ORFs (first and second frame) but one internal in the first
        alignment_dict = {
            'seq1': 'ATGAAAATGAAATAACCCTATGGGGTAG',
            'seq2': 'ATGAAAATGAAATAACCCTATGGGGTAG',
            'seq3': 'ATGAAAATGAAATAACCCTATGGGGTAG'
        }
        aln_fw = MSA(create_alignment(alignment_dict))
        aln_rw = MSA(create_alignment(aln_fw.calc_reverse_complement_alignment()))

        for aln in [aln_rw, aln_fw]:
            orfs = aln.get_non_overlapping_conserved_orfs(min_length=9)
            assert len(orfs) == 2

    def test_adjecent_orfs(self):
        """Test if ORFs that are directly adjacent (touching boundaries) are retained as non-overlapping for both frames."""
        alignment_dict = {
            'seq1': 'ATGAAAATGTAAATGAAAATGTAA',
            'seq2': 'ATGAAAATGTAAATGAAAATGTAA',
            'seq3': 'ATGAAAATGTAAATGAAAATGTAA'
        }
        aln_fw = MSA(create_alignment(alignment_dict))
        aln_rw = MSA(create_alignment(aln_fw.calc_reverse_complement_alignment()))
        for aln in [aln_rw, aln_fw]:
            orfs = aln.get_non_overlapping_conserved_orfs(min_length=9)
            assert len(orfs) == 2

    def test_non_overlapping_preserves_structure(self):
        """Test that non-overlapping ORFs preserve the same data structure."""
        alignment_dict = {
            'seq1': 'ATGAAATAA',
            'seq2': 'ATGAAATAA',
            'seq3': 'ATGAAATAA'
        }
        aln = MSA(create_alignment(alignment_dict))
        orfs = aln.get_non_overlapping_conserved_orfs(min_length=9)
        # first orf
        orf = orfs[list(orfs.keys())[0]]

        # Check all required keys are present
        assert 'location' in orf
        assert 'frame' in orf
        assert 'strand' in orf
        assert 'conservation' in orf
        assert 'internal' in orf

    def test_non_overlapping_with_identity_cutoff(self):
        """Test non-overlapping ORFs with identity cutoff."""
        alignment_dict = {
            'seq1': 'ATGAAATAAATGCCCTAG',
            'seq2': 'ATGAAATAAATGCCCTAG',
            'seq3': 'ATGAAATAAATGCCCTAG'
        }
        aln = MSA(create_alignment(alignment_dict))
        orfs = aln.get_non_overlapping_conserved_orfs(min_length=9, identity_cutoff=90.0)

        # All returned ORFs should meet identity cutoff
        for orf in orfs.values():
            assert orf['conservation'] >= 90.0
