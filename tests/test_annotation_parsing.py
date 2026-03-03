from pathlib import Path

from Bio import SeqIO
from msaexplorer.explore import Annotation, MSA


def test_annotation_parses_genbank_from_file() -> None:
    alignment = MSA("example_alignments/DNA.fasta")
    annotation = Annotation(alignment, "example_alignments/DNA_RNA.gb")

    assert annotation.ann_type == "gb"
    assert annotation.features
    assert "source" in annotation.features

    # Ensure expected feature shape is retained after parsing and mapping.
    first_source = annotation.features["source"][0]
    assert "location" in first_source
    assert "strand" in first_source
    assert isinstance(first_source["location"], list)


def test_annotation_parses_genbank_from_string() -> None:
    alignment = MSA("example_alignments/DNA.fasta")
    gb_content = Path("example_alignments/DNA_RNA.gb").read_text(encoding="utf-8")
    annotation = Annotation(alignment, gb_content)

    assert annotation.ann_type == "gb"
    assert annotation.features

    cds_features = annotation.features.get("CDS", {})
    assert cds_features
    first_cds = next(iter(cds_features.values()))
    assert "product" in first_cds


def test_annotation_parses_genbank_from_seqio_iterator() -> None:
    """Test that Annotation can accept a Bio.SeqIO iterator directly."""
    alignment = MSA("example_alignments/DNA.fasta")

    # Create a SeqIO iterator
    seqio_iterator = SeqIO.parse("example_alignments/DNA_RNA.gb", "genbank")

    # Pass the iterator directly to Annotation
    annotation = Annotation(alignment, seqio_iterator)

    assert annotation.ann_type == "gb"
    assert annotation.features
    assert "source" in annotation.features

    # Verify features are properly parsed and mapped
    cds_features = annotation.features.get("CDS", {})
    assert cds_features
    first_cds = next(iter(cds_features.values()))
    assert "product" in first_cds
    assert "location" in first_cds
    assert isinstance(first_cds["location"], list)


def test_annotation_parses_bed_from_file() -> None:
    """Test that Annotation can read BED format files."""
    alignment = MSA("example_alignments/DNA.fasta")
    annotation = Annotation(alignment, "tests/data/annotation.bed")

    assert annotation.ann_type == "bed"
    assert annotation.features
    assert "region" in annotation.features

    # Verify feature structure
    region_features = annotation.features["region"]
    assert len(region_features) >= 1
    first_region = region_features[0]
    assert "location" in first_region
    assert "strand" in first_region
    assert isinstance(first_region["location"], list)


def test_annotation_parses_bed_from_string() -> None:
    """Test that Annotation can read BED format from raw string."""
    alignment = MSA("example_alignments/DNA.fasta")
    bed_content = Path("tests/data/annotation.bed").read_text(encoding="utf-8")
    annotation = Annotation(alignment, bed_content)

    assert annotation.ann_type == "bed"
    assert annotation.features
    assert "region" in annotation.features

    region_features = annotation.features["region"]
    assert len(region_features) >= 1


def test_annotation_parses_gff3_from_file() -> None:
    """Test that Annotation can read GFF3 format files."""
    alignment = MSA("example_alignments/DNA.fasta")
    annotation = Annotation(alignment, "tests/data/annotation.gff3")

    assert annotation.ann_type == "gff"
    assert annotation.features
    # GFF3 parser converts 'region' to 'source' to match GenBank conventions
    assert "source" in annotation.features or "region" in annotation.features

    # Verify feature structure
    features_dict = annotation.features.get("source", {}) or annotation.features.get("region", {})
    assert len(features_dict) >= 1
    first_feature = next(iter(features_dict.values()))
    assert "location" in first_feature
    assert "strand" in first_feature


def test_annotation_parses_gff3_from_string() -> None:
    """Test that Annotation can read GFF3 format from raw string."""
    alignment = MSA("example_alignments/DNA.fasta")
    gff_content = Path("tests/data/annotation.gff3").read_text(encoding="utf-8")
    annotation = Annotation(alignment, gff_content)

    assert annotation.ann_type == "gff"
    assert annotation.features

    # GFF3 parser converts 'region' to 'source' to match GenBank conventions
    features_dict = annotation.features.get("source", {}) or annotation.features.get("region", {})
    assert len(features_dict) >= 1


