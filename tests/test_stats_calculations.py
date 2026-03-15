"""Tests for alignment statistics calculation methods in ``MSA``."""

import pytest
import numpy as np
from conftest import create_alignment
from msaexplorer.explore import MSA
from msaexplorer._msa_data_classes import PairwiseDistanceResult


class TestCalcEntropy:
    """Tests for calc_entropy."""

    def test_returns_list_of_correct_length(self):
        """Entropy list has length equal to alignment length."""
        msa = MSA(create_alignment({"s1": "ACGTACGT", "s2": "ACGTACGT"}))
        entropy = msa.calc_entropy()

        assert isinstance(entropy, list)
        assert len(entropy) == 8

    def test_identical_sequences_is_zero(self):
        """Identical sequences have zero entropy (no variation)."""
        msa = MSA(create_alignment({"s1": "ACGTACGT", "s2": "ACGTACGT"}))
        entropy = msa.calc_entropy()

        assert all(e == 0 for e in entropy)

    def test_perfectly_mixed_position(self):
        """Position with equal frequencies has maximum entropy."""
        msa = MSA(create_alignment({"s1": "AAAA", "s2": "CCCC", "s3": "GGGG", "s4": "TTTT"}))
        entropy = msa.calc_entropy()

        assert all(e == 1 for e in entropy)

    def test_normalized_between_zero_and_one(self):
        """Entropy values are normalized to [0, 1]."""
        msa = MSA(create_alignment({"s1": "ACGTACGT", "s2": "ACGTACGT", "s3": "GCGTACGT"}))
        entropy = msa.calc_entropy()

        assert all(0 <= e <= 1 for e in entropy)

    def test_single_gap_position(self):
        """Gaps are handled correctly (ignored in entropy calculation)."""
        msa = MSA(create_alignment({"s1": "A-GT", "s2": "ACGT", "s3": "ACGT"}))
        entropy = msa.calc_entropy()

        assert entropy[1] == 0

    def test_with_ambiguous_nucleotides(self):
        """Ambiguous nucleotides are handled in entropy calculation."""
        msa = MSA(create_alignment({"s1": "ARGT", "s2": "ACGT"}))
        entropy = msa.calc_entropy()

        assert entropy[1] == 0.75


class TestCalcGC:
    """Tests for calc_gc."""

    def test_returns_list_of_correct_length(self):
        """GC list has length equal to alignment length."""
        msa = MSA(create_alignment({"s1": "ACGTACGT", "s2": "ACGTACGT"}))
        gc = msa.calc_gc()

        assert isinstance(gc, list)
        assert len(gc) == 8

    def test_mixed_bases_partial_gc(self):
        """Mixed bases contribute to partial GC content."""
        msa = MSA(create_alignment({"s1": "AAGT", "s2": "ACGT"}))
        gc = msa.calc_gc()

        assert gc[0] == 0.0
        assert gc[1] == 0.5
        assert gc[2] == 1.0
        assert gc[3] == 0.0

    def test_raises_error_for_aa_alignment(self):
        """GC calculation raises TypeError for amino acid alignments."""
        msa = MSA(create_alignment({"s1": "ARNDCQ", "s2": "ARNDCQ"}))
        with pytest.raises(TypeError):
            msa.calc_gc()

    def test_with_ambiguous_nucleotides(self):
        """Ambiguous nucleotides contribute fractionally to GC."""
        msa = MSA(create_alignment({"s1": "SWGT", "s2": "ACGT"}))
        gc = msa.calc_gc()

        assert all(g == 0.5 for g in gc[0:1])


class TestCalcCoverage:
    """Tests for calc_coverage."""

    def test_returns_list_of_correct_length(self):
        """Coverage list has length equal to alignment length."""
        msa = MSA(create_alignment({"s1": "ACGTACGT", "s2": "ACGTACGT"}))
        coverage = msa.calc_coverage()

        assert isinstance(coverage, list)
        assert len(coverage) == 8

    def test_with_gaps(self):
        """Positions with all gaps have 0% coverage."""
        msa = MSA(create_alignment({"s1": "AT--", "s2": "C---", "s3": "G---", "s4": "GT--"}))
        coverage = msa.calc_coverage()

        assert coverage[0] == 1.0
        assert coverage[1] == 0.5
        assert all(c == 0.0 for c in coverage[2:])


class TestCalcTransitionTransversionScore:
    """Tests for calc_transition_transversion_score."""

    def test_returns_list_of_correct_length(self):
        """TS/TV score list has length equal to alignment length."""
        msa = MSA(create_alignment({"s1": "ACGTACGT", "s2": "ACGTACGT"}))
        score = msa.calc_transition_transversion_score()

        assert isinstance(score, list)
        assert len(score) == 8

    def test_identical_sequences_is_zero(self):
        """Identical sequences have zero TS/TV score."""
        msa = MSA(create_alignment({"s1": "ACGTACGT", "s2": "ACGTACGT"}))
        score = msa.calc_transition_transversion_score()

        assert all(s == 0 for s in score)

    def test_transition_is_positive(self):
        """Transition substitutions (A<->G, C<->T) are positive."""
        msa = MSA(create_alignment({"s1": "AGAG", "s2": "AAAA"}))
        score = msa.calc_transition_transversion_score()

        assert score[1] == 0.5
        assert score[3] == 0.5

    def test_transversion_is_negative(self):
        """Transversion substitutions (A<->C, A<->T, G<->C, G<->T) are negative."""
        msa = MSA(create_alignment({"s1": "ACAC", "s2": "AAAA"}))
        score = msa.calc_transition_transversion_score()

        assert score[1] == -0.5
        assert score[3] == -0.5


    def test_raises_error_for_aa_alignment(self):
        """TS/TV calculation raises TypeError for amino acid alignments."""
        msa = MSA(create_alignment({"s1": "ARNDCQ", "s2": "ARNDCQ"}))
        with pytest.raises(TypeError):
            msa.calc_transition_transversion_score()

    def test_rna_alignment(self):
        """TS/TV works with RNA alignments."""
        msa = MSA(create_alignment({"s1": "AGAGAG", "s2": "AUAAAA"}))
        score = msa.calc_transition_transversion_score()

        assert score[1] == -0.5


class TestGetSnps:
    """Tests for get_snps."""

    def test_no_variants(self):
        """Test if no variants are present"""
        msa = MSA(create_alignment({"s1": "ACGT", "s2": "ACGT"}))
        snps = msa.get_snps()

        assert snps["#CHROM"] == "consensus"
        assert snps["POS"] == {}

    def test_exclude_ambiguous(self):
        """Test if no variants are present for include_ambig=False"""
        msa = MSA(
            create_alignment({
                "ref": "AAAA",
                "q1": "ANAA",  # only ambiguous change at position 1
                "q2": "AAAA",
            }),
            reference_id="ref",
        )
        snps = msa.get_snps(include_ambig=False)

        assert snps["POS"] == {}

    def test_include_ambiguous(self):
        """Test if variants are present in correct frequency for include_ambig=False"""
        msa = MSA(
            create_alignment({
                "ref": "AAAA",
                "q1": "ACAA",  # non-ambiguous ALT C at position 1
                "q2": "ANAA",  # ambiguous ALT N at position 1
            }),
            reference_id="ref",
        )

        snps = msa.get_snps(include_ambig=True)
        alts = snps["POS"][1]["ALT"]

        assert 1 in snps["POS"]
        assert set(alts.keys()) == {"C", "N"}
        assert alts["C"]["AF"] == 0.5
        assert alts["N"]["AF"] == 0.5

    def test_alt_gaps_exclude_ambiguous(self):
        """Test that no gaps are present for include_ambig=False"""
        msa = MSA(
            create_alignment({
                "ref": "AAAA",
                "q1": "A-AA",  # only gap ALT at position 1
                "q2": "AAAA",
            }),
            reference_id="ref",
        )

        snps = msa.get_snps(include_ambig=False)

        assert snps["POS"] == {}

    def test_alt_gaps_include_ambiguous(self):
        """Test that there are gaps present for include_ambig=True"""
        msa = MSA(
            create_alignment({
                "ref": "AAAA",
                "q1": "A-AA",  # only gap ALT at position 1
                "q2": "ACAA",
            }),
            reference_id="ref",
        )

        snps = msa.get_snps(include_ambig=True)
        alts = snps["POS"][1]["ALT"]

        assert 1 in snps["POS"]
        assert set(alts.keys()) == {'-', 'C'}
        assert alts["-"]["AF"] == 0.5

    def test_gaps_in_reference(self):
        """Test that gaps are always included if in the reference"""
        msa = MSA(
            create_alignment({
                "ref": "A--A",
                "q1": "ACAA",  # non-ambiguous ALT C at position 1
                "q2": "A-AA",  # ambiguous ALT - at position 1
            }),
            reference_id="ref",
        )

        snps = msa.get_snps(include_ambig=False)

        assert snps["POS"][1]['ref'] == '-'
        assert snps["POS"][2]['ref'] == '-'

    def test_correct_sequence_identifiers(self):
        msa = MSA(
            create_alignment({
                "ref": "AAAA",
                "q1": "ACAA",
                "q2": "AGAA",
            }),
            reference_id="ref",
        )
        snps = msa.get_snps()

        assert snps["#CHROM"] == "ref"
        assert set(snps["POS"][1]["ALT"]["C"]["SEQ_ID"]) == {"q1"}
        assert set(snps["POS"][1]["ALT"]["G"]["SEQ_ID"]) == {"q2"}


class TestCalcPairwiseIdentityMatrix:
    """Tests for calc_pairwise_identity_matrix."""

    @pytest.mark.parametrize(
        "distance_type, expected_offdiag",
        [
            ("ghd", 57.14285714285714),
            ("lhd", 75.0),
            ("ged", 100.0),
            ("gcd", 50.0),
        ],
    )
    def test_distance_types(self, distance_type, expected_offdiag):
        msa = MSA(create_alignment({"s1": "-A-CGT-", "s2": "TAACG--"}))
        matrix = msa.calc_pairwise_identity_matrix(distance_type=distance_type)

        assert matrix.shape == (2, 2)
        assert matrix[0, 0] == 100.0
        assert matrix[1, 1] == 100.0
        assert matrix[0, 1] == pytest.approx(expected_offdiag)
        assert matrix[1, 0] == pytest.approx(expected_offdiag)

    def test_invalid_distance_type_raises(self):
        msa = MSA(create_alignment({"s1": "ACGT", "s2": "ACGT"}))
        with pytest.raises(ValueError, match="Invalid distance type"):
            msa.calc_pairwise_identity_matrix(distance_type="invalid")


class TestCalcPairwiseDistanceToReference:
    """Tests for calc_pairwise_distance_to_reference."""

    @pytest.mark.parametrize(
        "sequences",
        [
            {"ref": "AAAA", "q1": "AAAT", "q2": "AATT"},
            {"q1": "AAAT", "ref": "AAAA", "q2": "AATT"},
            {"q1": "AAAT", "q2": "AATT", "ref": "AAAA"},
        ],
    )
    def test_returns_dataclass_for_all_reference_positions(self, sequences):
        msa = MSA(create_alignment(sequences), reference_id="ref")

        result = msa.calc_pairwise_distance_to_reference(distance_type="ghd")

        assert isinstance(result, PairwiseDistanceResult)
        assert result.reference_id == "ref"
        assert result.sequence_ids == ["q1", "q2"]
        assert np.allclose(result.distances, np.array([75.0, 50.0]))

    def test_uses_consensus_when_reference_id_is_not_set(self):
        msa = MSA(create_alignment({"q1": "AAAA", "q2": "AAAT", "q3": "AATT"}))

        result = msa.calc_pairwise_distance_to_reference(distance_type="ghd")

        assert result.reference_id == "consensus"
        assert result.sequence_ids == ["q1", "q2", "q3"]
        assert np.allclose(result.distances, np.array([75.0, 100.0, 75.0]))

    def test_invalid_distance_type_raises(self):
        msa = MSA(create_alignment({"s1": "ACGT", "s2": "ACGT"}))

        with pytest.raises(ValueError, match="Invalid distance type"):
            msa.calc_pairwise_distance_to_reference(distance_type="invalid")


class TestCalcPositionMatrix:
    """Tests for calc_position_matrix."""

    expected = np.array([[2, 0, 0, 0], # A
                         [0, 0, 1, 2], # T
                         [0, 0, 1, 0], # G
                         [0, 2, 0, 0]]) # C

    def test_pfm_has_expected_counts(self):
        """Test if frequency counts in PFM are correct."""
        msa = MSA(create_alignment({"s1": "ACTT", "s2": "ACGT"}))
        pfm = msa.calc_position_matrix("PFM")

        assert (pfm == self.expected).all()

    def test_invalid_matrix_type_raises(self):
        """Test that invalid matrix type raises ValueError."""
        msa = MSA(create_alignment({"s1": "ACGT", "s2": "ACGT"}))
        with pytest.raises(ValueError, match="Matrix_type must be PFM, PPM, IC or PWM"):
            msa.calc_position_matrix("INVALID")

