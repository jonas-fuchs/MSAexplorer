"""Tests for the main plotting helpers in ``msaexplorer.draw``."""

from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pytest
from matplotlib.axes import Axes

from conftest import create_alignment
from msaexplorer import draw
from msaexplorer.explore import MSA


@pytest.fixture
def ax() -> Axes:
    """Create a fresh matplotlib axis for each test."""
    fig, axis = plt.subplots()
    yield axis
    plt.close(fig)


@pytest.fixture
def dna_msa() -> MSA:
    """Create a small nucleotide alignment with mismatches, a gap and a mask character."""
    return MSA(
        create_alignment(
            {
                'ref': 'ATGAAATTTTAA',
                'seq1': 'ATGAAATTTTAA',
                'seq2': 'ATGACATTTTAA',
                'seq3': 'ATGAAAT-TNAA',
            }
        ),
        reference_id='ref',
    )


@pytest.fixture
def aa_msa() -> MSA:
    """Create a small amino-acid alignment for AA plotting tests."""
    return MSA(
        create_alignment(
            {
                'ref_aa': 'MKTWQALV',
                'seq1_aa': 'MKTWQALV',
                'seq2_aa': 'MKAWQALV',
                'seq3_aa': 'MKTWQ-LV',
            }
        ),
        reference_id='ref_aa',
    )


@pytest.fixture
def orf_msa() -> MSA:
    """Create a small nucleotide alignment with one conserved ORF."""
    return MSA(
        create_alignment(
            {
                'ref': 'ATGAAATTTTAA',
                'seq1': 'ATGAAATTTTAA',
                'seq2': 'ATGAAGTTTTAA',
            }
        ),
        reference_id='ref',
    )


@pytest.fixture
def bed_annotation_file(tmp_path: Path) -> Path:
    """Create a BED annotation file that matches the test alignment ids."""
    annotation_file = tmp_path / 'regions.bed'
    annotation_file.write_text(
        'ref\t0\t4\n'
        'ref\t6\t10\n',
        encoding='utf-8',
    )
    return annotation_file


def _assert_axis_is_populated(ax: Axes) -> None:
    """Assert that plotting added visible artists to the axis."""
    assert ax.collections or ax.patches or ax.lines or ax.texts


def test_alignment_returns_axis_and_adds_collections(dna_msa: MSA, ax: Axes) -> None:
    """Test that ``alignment`` plots onto the provided axis."""
    returned_ax = draw.alignment(dna_msa, ax=ax, show_sequence_all=False)

    assert returned_ax is ax
    assert len(ax.collections) >= 1
    _assert_axis_is_populated(ax)


def test_identity_alignment_returns_axis_and_adds_collections(dna_msa: MSA, ax: Axes) -> None:
    """Test that ``identity_alignment`` plots onto the provided axis."""
    returned_ax = draw.identity_alignment(dna_msa, ax=ax, show_identity_sequence=True)

    assert returned_ax is ax
    assert len(ax.collections) >= 1
    _assert_axis_is_populated(ax)


def test_similarity_alignment_returns_axis_and_adds_collections(aa_msa: MSA, ax: Axes) -> None:
    """Test that ``similarity_alignment`` plots onto the provided axis."""
    returned_ax = draw.similarity_alignment(aa_msa, ax=ax, show_similarity_sequence=True)

    assert returned_ax is ax
    assert len(ax.collections) >= 1
    _assert_axis_is_populated(ax)


@pytest.mark.parametrize('stat_type', ['entropy', 'identity'])
def test_stat_plot_returns_axis_for_main_stat_types(dna_msa: MSA, ax: Axes, stat_type: str) -> None:
    """Test that ``stat_plot`` renders both scalar and matrix-based statistics."""
    returned_ax = draw.stat_plot(dna_msa, stat_type=stat_type, ax=ax, rolling_average=1)

    assert returned_ax is ax
    assert len(ax.collections) >= 1
    _assert_axis_is_populated(ax)


def test_variant_plot_returns_axis_and_draws_variant_markers(dna_msa: MSA, ax: Axes) -> None:
    """Test that ``variant_plot`` draws lollipop markers for detected variants."""
    returned_ax = draw.variant_plot(dna_msa, ax=ax, show_legend=True)

    assert returned_ax is ax
    assert len(ax.lines) >= 1
    assert ax.get_ylabel() == 'reference'
    _assert_axis_is_populated(ax)


def test_annotation_plot_returns_axis_and_draws_feature_patches(dna_msa: MSA, ax: Axes, bed_annotation_file: Path) -> None:
    """Test that ``annotation_plot`` draws annotation boxes from a BED file."""
    returned_ax = draw.annotation_plot(dna_msa, str(bed_annotation_file), feature_to_plot='ignored', ax=ax)

    assert returned_ax is ax
    assert len(ax.patches) >= 1
    assert 'bed regions' in ax.get_title(loc='left')
    _assert_axis_is_populated(ax)


def test_orf_plot_returns_axis_and_draws_orf_patches(orf_msa: MSA, ax: Axes) -> None:
    """Test that ``orf_plot`` draws at least one conserved ORF."""
    returned_ax = draw.orf_plot(orf_msa, ax=ax, min_length=9, non_overlapping_orfs=True)

    assert returned_ax is ax
    assert len(ax.patches) >= 1
    assert ax.get_title(loc='left') == 'conserved orfs'
    _assert_axis_is_populated(ax)


@pytest.mark.parametrize('plot_type', ['stacked', 'logo'])
def test_sequence_logo_returns_axis_for_both_plot_modes(dna_msa: MSA, ax: Axes, plot_type: str) -> None:
    """Test that ``sequence_logo`` supports both stacked and logo rendering."""
    returned_ax = draw.sequence_logo(dna_msa, ax=ax, plot_type=plot_type)

    assert returned_ax is ax
    _assert_axis_is_populated(ax)


def test_consensus_plot_returns_axis_and_draws_collection_and_text(dna_msa: MSA, ax: Axes) -> None:
    """Test that ``consensus_plot`` draws a consensus row on the provided axis."""
    returned_ax = draw.consensus_plot(dna_msa, ax=ax, show_sequence=True)

    assert returned_ax is ax
    assert len(ax.collections) >= 1
    assert len(ax.texts) == dna_msa.length
    _assert_axis_is_populated(ax)


def test_simplot_returns_axis_and_draws_one_line_per_non_reference_sequence(dna_msa: MSA, ax: Axes) -> None:
    """Test that ``simplot`` draws one trace per non-reference sequence."""
    returned_ax = draw.simplot(
        dna_msa,
        ref='ref',
        ax=ax,
        window_size=4,
        step_size=2,
        distance_calculation='k2p',
        show_legend=True,
        show_reference=True,
    )

    assert returned_ax is ax
    assert len(ax.lines) == len(dna_msa) - 1
    assert ax.get_ylabel() == 'similarity (%)'
    _assert_axis_is_populated(ax)

