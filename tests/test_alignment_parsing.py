import pathlib
import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment

from msaexplorer.explore import MSA
from msaexplorer._helpers import _read_alignment

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
    parsed_alignment = _read_alignment(str(alignment_file))

    assert parsed_alignment == EXPECTED_ALIGNMENT


@pytest.mark.parametrize("format_name", ["fasta", "clustal", "phylip", "nexus", "stockholm"])
def test_read_alignment_parses_multiple_formats_from_string(format_name: str):
    """Test that MSA can read different alignment formats from a raw string."""
    alignment_file = DATA_DIR / f"alignment.{format_name}"
    alignment_content = alignment_file.read_text(encoding="utf-8")
    parsed_alignment = _read_alignment(alignment_content)

    assert parsed_alignment == EXPECTED_ALIGNMENT


def test_read_alignment_accepts_bio_alignment_object():
    """Test that MSA can accept Bio.Align.MultipleSeqAlignment objects directly."""
    # Create a Bio alignment object
    records = [
        SeqRecord(Seq(sequence), id=seq_id, description="")
        for seq_id, sequence in EXPECTED_ALIGNMENT.items()
    ]
    bio_alignment = MultipleSeqAlignment(records)

    # Pass it directly to MSA._read_alignment
    parsed_alignment = _read_alignment(bio_alignment)

    assert parsed_alignment == EXPECTED_ALIGNMENT


def test_msa_initializes_with_bio_alignment_object():
    """Test that MSA can be initialized with a Bio.Align.MultipleSeqAlignment object."""
    # Create a Bio alignment object
    records = [
        SeqRecord(Seq(sequence), id=seq_id, description="")
        for seq_id, sequence in EXPECTED_ALIGNMENT.items()
    ]
    bio_alignment = MultipleSeqAlignment(records)

    # Initialize MSA with Bio alignment object
    msa = MSA(bio_alignment)

    assert msa.alignment == EXPECTED_ALIGNMENT
    assert msa.aln_type in ["DNA", "RNA", "AA"]
    assert len(msa.alignment) == 3


def test_read_alignment_raises_for_unparseable_content() -> None:
    """Test that MSA raises ValueError when given unparseable content."""
    with pytest.raises(ValueError, match="could not be parsed"):
        _read_alignment("this is not an alignment")


def test_zoom_allows_exclusive_end_equal_to_alignment_length() -> None:
    """Accept zoom ranges that use Python slicing semantics [start:end)."""
    msa = MSA(str(DATA_DIR / 'alignment.fasta'))

    msa.zoom = (0, msa.length)

    assert msa.length == len(EXPECTED_ALIGNMENT['seq1'])


def test_zoom_rejects_end_beyond_alignment_length() -> None:
    """Reject zoom ranges whose end exceeds the alignment length."""
    msa = MSA(str(DATA_DIR / 'alignment.fasta'))

    with pytest.raises(ValueError, match='Zoom position out of range'):
        msa.zoom = (0, msa.length + 1)

