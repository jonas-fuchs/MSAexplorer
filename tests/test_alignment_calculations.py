"""Tests for alignment calculation methods in ``MSA``."""

import numpy as np
import pytest

from conftest import create_alignment
from msaexplorer.explore import MSA


class TestCalcNumericalAlignment:
    def test_returns_expected_shape(self):
        msa = MSA(create_alignment({"s1": "ACGT", "s2": "ACGT"}))
        matrix = msa.calc_numerical_alignment()
        assert matrix.shape == (2, 4)

    def test_identical_sequences_have_identical_rows(self):
        msa = MSA(create_alignment({"s1": "ACGT", "s2": "ACGT"}))
        matrix = msa.calc_numerical_alignment()
        assert np.array_equal(matrix[0], matrix[1])

    def test_encode_mask_marks_n_as_minus_two(self):
        msa = MSA(create_alignment({"s1": "ANGT", "s2": "ACGT"}))
        matrix = msa.calc_numerical_alignment(encode_mask=True)
        assert matrix[0, 1] == -2

    def test_encode_ambiguities_marks_iupac_as_minus_three(self):
        msa = MSA(create_alignment({"s1": "ARGT", "s2": "ACGT"}))
        matrix = msa.calc_numerical_alignment(encode_ambiguities=True)
        assert matrix[0, 1] == -3


class TestCalcIdentityAlignment:
    def test_basic_identity_and_mismatch_encoding(self):
        msa = MSA(create_alignment({"s1": "ACGT", "s2": "ACGT", "s3": "ACCT"}))
        matrix = msa.calc_identity_alignment()

        assert matrix.shape == (3, 4)
        assert matrix[0, 0] == 0
        assert matrix[2, 2] == -1

    def test_gap_is_nan_by_default(self):
        msa = MSA(create_alignment({"s1": "AC-T", "s2": "ACGT"}))
        matrix = msa.calc_identity_alignment()
        assert np.isnan(matrix[0, 2])

    def test_encode_gaps_false_keeps_gap_non_nan(self):
        msa = MSA(create_alignment({"s1": "AC-T", "s2": "ACGT"}))
        matrix = msa.calc_identity_alignment(encode_gaps=False)
        assert not np.isnan(matrix[0, 2])

    def test_reference_id_is_used(self):
        msa = MSA(create_alignment({"ref": "AAAA", "q1": "AAAT", "q2": "AATA"}), reference_id="ref")
        matrix = msa.calc_identity_alignment()
        assert matrix[1, 3] == -1
        assert matrix[2, 2] == -1
        assert np.all(matrix[0] == 0)

    def test_encode_each_mismatch_char_avoids_minus_one(self):
        msa = MSA(create_alignment({"s1": "AAAA", "s2": "AAAT"}))
        matrix = msa.calc_identity_alignment(encode_each_mismatch_char=True)
        vals = matrix[~np.isnan(matrix)]
        assert -1 not in vals
        assert matrix[1, 3] > 0


class TestCalcSimilarityAlignment:
    def test_normalized_values_in_unit_interval(self):
        msa = MSA(create_alignment({"s1": "ACTA", "s2": "ACCA", "s3": "ACTT"}))
        sim = msa.calc_similarity_alignment(normalize=True)
        vals = sim[~np.isnan(sim)]
        assert np.all(vals >= 0)
        assert np.all(vals <= 1)

    def test_default_matrix_for_dna_is_trans(self):
        msa = MSA(create_alignment({"s1": "ACGT", "s2": "ACGT"}))
        sim = msa.calc_similarity_alignment()
        assert sim.dtype.metadata["matrix"] == "TRANS"

    def test_default_matrix_for_aa_is_blosum65(self):
        msa = MSA(create_alignment({"s1": "ARNDCQ", "s2": "ARNDCQ"}))
        sim = msa.calc_similarity_alignment()
        assert sim.dtype.metadata["matrix"] == "BLOSUM65"

    def test_invalid_matrix_raises_value_error(self):
        msa = MSA(create_alignment({"s1": "ACGT", "s2": "ACGT"}))
        with pytest.raises(ValueError):
            msa.calc_similarity_alignment(matrix_type="NOPE")

    def test_reference_gap_column_gets_score_one_for_non_reference_rows(self):
        msa = MSA(create_alignment({"ref": "A-TA", "q1": "ACTA", "q2": "AGTA"}), reference_id="ref")
        sim = msa.calc_similarity_alignment(normalize=True)

        assert np.isnan(sim[0, 1])
        assert np.isclose(sim[1, 1], 1.0)
        assert np.isclose(sim[2, 1], 1.0)

