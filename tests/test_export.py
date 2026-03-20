"""Tests for export helpers in ``msaexplorer.export``."""

import numpy as np
import pytest

from msaexplorer import export
from msaexplorer._data_classes import AlignmentStats, OpenReadingFrame, OrfCollection, SingleNucleotidePolymorphism, VariantCollection


@pytest.fixture
def snp_result():
    """Returns a SNP dataclass object for testing."""
    return VariantCollection(
        chrom="ref",
        positions={
            3: SingleNucleotidePolymorphism(
                ref="C",
                alt={"T": (1.0, ("s1",))},
            ),
            1: SingleNucleotidePolymorphism(
                ref="A",
                alt={
                    "G": (0.5, ("s1", "s2")),
                    "-": (0.5, ("s3",)),
                },
            ),
        },
    )


class TestSnpsExport:
    def test_vcf_output_is_correct_and_sorted(self, snp_result):
        """Test that the VCF has correct vcf format."""
        result = export.snps(snp_result, format_type="vcf")

        assert result == "\n".join(
            [
                "##fileformat=VCFv4.2",
                "##source=MSAexplorer",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
                "ref\t2\t.\tA\tG,-\t.\t.\tAF=0.5,0.5;SEQ_ID=s1|s2,s3",
                "ref\t4\t.\tC\tT\t.\t.\tAF=1.0;SEQ_ID=s1",
            ]
        )

    def test_tabular_output_is_correct(self, snp_result):
        """Test that the VCF has correct tabular format."""
        result = export.snps(snp_result, format_type="tabular")

        assert result == "\n".join(
            [
                "CHROM\tPOS\tREF\tALT\tAF\tSEQ_ID",
                "ref\t2\tA\tG\t0.5\ts1,s2",
                "ref\t2\tA\t-\t0.5\ts3",
                "ref\t4\tC\tT\t1.0\ts1",
            ]
        )

    def test_invalid_input_raises_value_error(self):
        """Test that invalid input raises a ValueError."""
        with pytest.raises(ValueError, match="must be a VariantCollection"):
            export.snps([], format_type="vcf")

    def test_invalid_format_raises_value_error(self, snp_result):
        with pytest.raises(ValueError, match="Invalid format_type"):
            export.snps(snp_result, format_type="csv")

    def test_file_export_creates_expected_extension(self, snp_result, tmp_path):
        """Test that the file extension is correct when written."""
        path_without_ext = tmp_path / "nested" / "snps_output"

        result = export.snps(snp_result, format_type="vcf", path=str(path_without_ext))

        assert result is None
        written = path_without_ext.with_suffix(path_without_ext.suffix + ".vcf")
        assert written.exists()
        assert "##fileformat=VCFv4.2" in written.read_text()


class TestFastaExport:
    def test_single_sequence_to_string(self):
        """Test that a single sequence is exported correctly."""
        result = export.fasta("ACGT", header="s1")

        assert result == ">s1\nACGT"

    def test_dictionary_sequence_to_string(self):
        """Test that a dictionary of sequences is exported correctly."""
        result = export.fasta({"s1": "ACGT", "s2": "A-GT"})

        assert result == ">s1\nACGT\n>s2\nA-GT"

    def test_invalid_sequence_raises(self):
        """Test that invalid sequences raise a ValueError."""
        with pytest.raises(ValueError, match="invalid characters"):
            export.fasta("ACGTZ", header="s1")

    def test_file_export_writes_content(self, tmp_path):
        """Test that the file is written correctly."""
        out = tmp_path / "subdir" / "seq.fasta"

        result = export.fasta("ACGT", header="s1", path=str(out))

        assert result is None
        assert out.exists()
        assert out.read_text() == ">s1\nACGT"


class TestStatsExport:
    def test_alignment_stats_dataclass_output(self):
        """Test that the AlignmentStats dataclass is exported correctly."""
        stat_data = AlignmentStats(
            stat_name="entropy",
            positions=np.array([2, 4, 6]),
            values=np.array([0.1, 0.2, 0.3]),
        )

        result = export.stats(stat_data, seperator=",")

        assert result == "\n".join([
            "position,value",
            "2,0.1",
            "4,0.2",
            "6,0.3",
        ])

    def test_plain_array_output_uses_zero_based_positions(self):
        """Test that plain arrays are exported correctly."""
        result = export.stats([10, 20, 30])

        assert result == "\n".join([
            "position\tvalue",
            "0\t10",
            "1\t20",
            "2\t30",
        ])


class TestOrfExport:
    @pytest.fixture
    def orf_collection(self):
        return OrfCollection(orfs=(
            OpenReadingFrame(
                orf_id='orf_1',
                location=((10, 50),),
                frame=1,
                strand='+',
                conservation=87.125,
            ),
        ))

    def test_orf_output_is_correct(self, orf_collection):
        result = export.orf(orf_collection, chrom="chr1")

        assert result == "chr1\t10\t50\torf_1\t87.12\t+"

    def test_orf_invalid_input_raises(self):
        with pytest.raises(ValueError, match="instance"):
            export.orf(OrfCollection(), chrom="chr1")
            export.orf({}, chrom="chr1")


class TestCharacterFreqExport:
    def test_character_freq_output_is_correct(self):
        """Test that the character frequency output is correct."""
        data = {
            "total": {},
            "seq1": {
                "A": {"counts": 2, "% of alignment": 50.0, "% of non-gapped": 66.67},
                "C": {"counts": 1, "% of alignment": 25.0, "% of non-gapped": 33.33},
                "-": {"counts": 1, "% of alignment": 25.0, "% of non-gapped": 0.0},
            },
        }

        result = export.character_freq(data, seperator=",")

        assert result == "\n".join([
            "sequence,char,counts,% of non-gapped",
            "seq1,A,2,66.67",
            "seq1,C,1,33.33",
        ])

    def test_character_freq_invalid_char_raises(self):
        """Test that invalid characters raise a ValueError."""
        bad = {"seq1": {"?": {"counts": 1, "% of alignment": 10.0, "% of non-gapped": 10.0}}}

        with pytest.raises(ValueError, match="invalid"):
            export.character_freq(bad)


class TestPercentRecoveryExport:
    def test_percent_recovery_output_is_correct(self):
        """Test that the percent recovery output is correct."""
        rec = {"seq1": 75.0, "seq2": 90.5}

        result = export.percent_recovery(rec)

        assert result == "\n".join([
            "sequence\t% recovery",
            "seq1\t75.0",
            "seq2\t90.5",
        ])

    def test_percent_recovery_value_must_be_float(self):
        """Test that the percent recovery value must be a float."""
        with pytest.raises(ValueError, match="invalid"):
            export.percent_recovery({"seq1": 75})
