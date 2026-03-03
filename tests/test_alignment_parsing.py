import pathlib
import pytest
from msaexplorer.explore import MSA

DATA_DIR = pathlib.Path(__file__).parent / "data"

EXPECTED_ALIGNMENT = {
    "seq1": "ATGCATGC",
    "seq2": "ATG-ATGC",
    "seq3": "ATGCAT-T",
}


@pytest.mark.parametrize("format_name", ["fasta", "clustal", "phylip", "nexus", "stockholm"])
def test_read_alignment_parses_multiple_formats_from_file(format_name: str):
    """Test that MSA can read different alignment formats from a file."""
    alignment_file = DATA_DIR / f"alignment.{format_name}"
    parsed_alignment = MSA._read_alignment(str(alignment_file))

    assert parsed_alignment == EXPECTED_ALIGNMENT


@pytest.mark.parametrize("format_name", ["fasta", "clustal", "phylip", "nexus", "stockholm"])
def test_read_alignment_parses_multiple_formats_from_string(format_name: str):
    """Test that MSA can read different alignment formats from a raw string."""
    alignment_file = DATA_DIR / f"alignment.{format_name}"
    alignment_content = alignment_file.read_text(encoding="utf-8")
    parsed_alignment = MSA._read_alignment(alignment_content)

    assert parsed_alignment == EXPECTED_ALIGNMENT


def test_read_alignment_raises_for_unparseable_content() -> None:
    """Test that MSA raises ValueError when given unparseable content."""
    with pytest.raises(ValueError, match="could not be parsed"):
        MSA._read_alignment("this is not an alignment")

