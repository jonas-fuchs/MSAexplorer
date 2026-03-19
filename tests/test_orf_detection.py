"""
Tests for conserved ORF detection.
"""

import pytest
from msaexplorer.explore import MSA
from msaexplorer._data_classes import OpenReadingFrame, OrfContainer
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
        orf_obj = orfs[0]
        assert orf_obj.strand == '+'
        assert orf_obj.location[0][0] == 0  # Start at position 0
        assert orf_obj.location[0][1] == 12  # End at position 12

    def test_with_gaps(self):
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
        """Test that empty collection is returned when no ORFs found."""
        alignment_dict = {
            'seq1': 'AAAAAAAAAAA',
            'seq2': 'AAAAAAAAAAA',
            'seq3': 'AAAAAAAAAAA'
        }
        aln = MSA(create_alignment(alignment_dict))
        orfs = aln.get_conserved_orfs(min_length=9)

        assert len(orfs) == 0

    def test_with_internal_stop(self):
        """Test that ORFs with internal stop codons are rejected."""
        alignment_dict = {
            'seq1': 'ATGTAAATGTAA',
            'seq2': 'ATGTAAATGTAA',
            'seq3': 'ATGTAAATGTAA'
        }
        aln = MSA(create_alignment(alignment_dict))
        orfs = aln.get_conserved_orfs(min_length=9)

        assert len(orfs) == 0  # Should not find ORF because of internal stop at position 3-5

    def test_multiple_frames(self):
        """Test ORF detection across different reading frames."""
        alignment_dict = {
            'seq1': 'CATGAAATAAGATGCCCTAG',
            'seq2': 'CATGAAATAAGATGCCCTAG',
            'seq3': 'CATGAAATAAGATGCCCTAG'
        }
        aln = MSA(create_alignment(alignment_dict))
        orfs = aln.get_conserved_orfs(min_length=9)

        # Should detect ORFs in frame 1
        frames = [orf.frame for orf in orfs.values()]
        assert 1 in frames

    def test_reverse_strand(self):
        """Test ORF detection on reverse complement strand."""
        alignment_dict = {
            'seq1': 'ATGTTATTTCATTAA',
            'seq2': 'ATGTTATTTCATTAA',
            'seq3': 'ATGTTATTTCATTAA'
        }
        aln = MSA(create_alignment(alignment_dict))
        orfs = aln.get_conserved_orfs(min_length=9)
        strands = [orf.strand for orf in orfs.values()]

        assert strands == ['+', '-']

    def test_non_conserved_start(self):
        """Test that non-conserved start codons are rejected."""
        alignment_dict = {
            'seq1': 'ATGAAATAA',
            'seq2': 'ATGAAATAA',
            'seq3': 'TTGAAATAA'   # Has TTG (not a start codon)
        }
        aln = MSA(create_alignment(alignment_dict))
        orfs = aln.get_conserved_orfs(min_length=9)

        assert len(orfs) == 0

    def test_non_conserved_stop(self):
        """Test that non-conserved stop codons are rejected."""
        alignment_dict = {
            'seq1': 'ATGAAATAA',
            'seq2': 'ATGAAATAA',
            'seq3': 'ATGAAACAA'   # Has CAA (not a stop)
        }
        aln = MSA(create_alignment(alignment_dict))
        orfs = aln.get_conserved_orfs(min_length=9)

        assert len(orfs) == 0

    def test_identity_cutoff(self):
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

    def test_internal_orfs_detected(self):
        """Test that internal ORFs are properly marked."""
        # Two start codons to same stop --> one internal
        alignment_dict = {
            'seq1': 'ATGAAAATGCCCTAA',
            'seq2': 'ATGAAAATGCCCTAA',
            'seq3': 'ATGAAAATGCCCTAA'
        }
        aln = MSA(create_alignment(alignment_dict))
        orfs = aln.get_conserved_orfs(min_length=9)

        assert len(orfs['ORF_0'].internal) == 1

    def test_rna_alignment(self):
        """Test ORF detection works with RNA alignments (U instead of T)."""
        alignment_dict = {
            'seq1': 'AUGAAAUAA',
            'seq2': 'AUGAAAUAA',
            'seq3': 'AUGAAAUAA'
        }
        aln = MSA(create_alignment(alignment_dict))
        orfs = aln.get_conserved_orfs(min_length=9)

        assert len(orfs) == 1

    def test_amino_acid_raises_error(self):
        """Test that amino acid alignments raise TypeError."""
        alignment_dict = {
            'seq1': 'ARNDCQEGHILKMFPSTWYV',
            'seq2': 'ARNDCQEGHILKMFPSTWYV',
            'seq3': 'ARNDCQEGHILKMFPSTWYV'
        }
        aln = MSA(create_alignment(alignment_dict))

        with pytest.raises(TypeError, match='ORF search only for RNA/DNA alignments'):
            aln.get_conserved_orfs(min_length=9)

    def test_invalid_min_length(self):
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

    def test_invalid_identity_cutoff(self):
        """Test that invalid identity_cutoff values raise ValueError."""
        alignment_dict = {
            'seq1': 'ATGAAATAA',
            'seq2': 'ATGAAATAA',
            'seq3': 'ATGAAATAA'
        }
        aln = MSA(create_alignment(alignment_dict))

        with pytest.raises(ValueError, match='conservation cutoff must be between 0 and 100'):
            aln.get_conserved_orfs(min_length=9, identity_cutoff=-1.0)

        with pytest.raises(ValueError, match='conservation cutoff must be between 0 and 100'):
            aln.get_conserved_orfs(min_length=9, identity_cutoff=101.0)

    def test_return_structure(self):
        """Test that the return value is an OrfContainer instance with correct structure."""
        alignment_dict = {
            'seq1': 'ATGAAATAA',
            'seq2': 'ATGAAATAA',
            'seq3': 'ATGAAATAA'
        }
        aln = MSA(create_alignment(alignment_dict))
        orfs = aln.get_conserved_orfs(min_length=9)

        assert isinstance(orfs, OrfContainer)

        orf_obj = orfs[0]
        assert isinstance(orf_obj, OpenReadingFrame)
        assert isinstance(orf_obj.location, tuple)
        assert isinstance(orf_obj.frame, int)
        assert orf_obj.strand in ['+', '-']
        assert isinstance(orf_obj.conservation, float)
        assert isinstance(orf_obj.internal, tuple)

    def test_getitem_by_name(self):
        """OrfContainer supports access by orf_id string."""
        alignment_dict = {
            'seq1': 'ATGAAATAA',
            'seq2': 'ATGAAATAA',
            'seq3': 'ATGAAATAA'
        }
        aln = MSA(create_alignment(alignment_dict))
        orfs = aln.get_conserved_orfs(min_length=9)

        assert orfs['ORF_0'] is orfs[0]

    def test_contains(self):
        """'ORF_0' should be contained in result, 'ORF_99' should not."""
        alignment_dict = {
            'seq1': 'ATGAAATAA',
            'seq2': 'ATGAAATAA',
            'seq3': 'ATGAAATAA'
        }
        aln = MSA(create_alignment(alignment_dict))
        orfs = aln.get_conserved_orfs(min_length=9)

        assert 'ORF_0' in orfs
        assert 'ORF_99' not in orfs

    def test_orf_len_and_contains_position(self):
        """OpenReadingFrame.__len__ and __contains__ work correctly."""
        alignment_dict = {
            'seq1': 'ATGAAATAA',
            'seq2': 'ATGAAATAA',
            'seq3': 'ATGAAATAA'
        }
        aln = MSA(create_alignment(alignment_dict))
        orfs = aln.get_conserved_orfs(min_length=9)
        orf_obj = orfs[0]

        assert len(orf_obj) == 9        # 0..9
        assert 0 in orf_obj             # start is inside
        assert 8 in orf_obj             # last position inside
        assert 9 not in orf_obj         # stop is exclusive


class TestGetNonOverlappingConservedOrfs:
    """Test the get_non_overlapping_conserved_orfs method."""

    def test_basic(self):
        """Test basic non-overlapping ORF selection."""
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
        """Test if ORFs that are directly adjacent are retained as non-overlapping."""
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

    def test_with_identity_cutoff(self):
        """Test non-overlapping ORFs with identity cutoff."""
        alignment_dict = {
            'seq1': 'ATGAAATAAATGCCCTAG',
            'seq2': 'ATGAAATAAATGCCCTAG',
            'seq3': 'ATGAAATAAATGCCCTAG'
        }
        aln = MSA(create_alignment(alignment_dict))
        orfs = aln.get_non_overlapping_conserved_orfs(min_length=9, identity_cutoff=90.0)

        for orf in orfs.values():
            assert orf.conservation >= 90.0
