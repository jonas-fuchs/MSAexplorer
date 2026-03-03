r"""
# MSAexplorer

MSAexplorer is a Python package for analyzing and visualizing multiple sequence alignments (MSAs).
It provides statistics, publication-quality plots, and flexible data export with minimal dependencies.

## Key Features

- **Multiple input formats**: FASTA, CLUSTAL, PHYLIP, STOCKHOLM, NEXUS, or Bio.Align.MultipleSeqAlignment
- **Annotation support**: GenBank, GFF3, BED with automatic coordinate mapping to alignment
- **Statistical analysis**: Entropy, GC content, coverage, pairwise identity, SNPs, ORF detection
- **Publication ready plots**: Identity matrices, statistical plots, sequence logos, variant plots, ORF tracks
- **Flexible export**: SNPs (VCF/tabular), FASTA, statistics, ORFs, recovery, character frequencies
- **Interactive web app**: Shiny interface (run locally or use online at GitHub Pages)
- **Biopython integration**: Seamless intergration with Biopython objects and formats

## Installation

```bash
pip install msaexplorer  # only the python package

# or with the optional app dependencies:
pip install msaexplorer[app]
pip install msaexplorer[app-plus]  # adds pyfamsa and pytrimal
```

## Quick Tutorial

### Basic Usage

```python
from msaexplorer import explore

# Load alignment (auto-detects format, supports FASTA, CLUSTAL, PHYLIP, STOCKHOLM, NEXUS)
msa = explore.MSA('alignment.fasta')
# or from Biopython object
from Bio import AlignIO
alignment = AlignIO.read('alignment.fasta', 'fasta')
msa = explore.MSA(alignment)

# Get alignment info
print(f"Sequences: {len(msa.alignment)}, Length: {msa.length}, Type: {msa.aln_type}")

# Compute statistics such as:
entropy = msa.calc_entropy()
gc = msa.calc_gc()
coverage = msa.calc_coverage()
```

### Load annotations and link to alignment

```python
from msaexplorer import explore

msa = explore.MSA('alignment.fasta')

# Load annotation (GenBank, GFF3, or BED - auto-detects format)
annotation = explore.Annotation(msa, 'annotation.gb')

# Or from Biopython iterator
from Bio import SeqIO
annotation = explore.Annotation(msa, SeqIO.parse('annotation.gb', 'genbank'))

# Access features
print(annotation.features.keys())  # Feature types
print(annotation.ann_type)          # 'gb', 'gff', or 'bed'
```

### Reference and Zoom

```python
msa = explore.MSA('alignment.fasta')

# Set reference for identity calculations
msa.reference_id = 'seq1'

# Focus on region
msa.zoom = (100, 500)  # or just msa.zoom = 100 for start at 0 to 100
```

### Analysis Functions (MSA class)

```python
# SNP analysis (DNA/RNA)
snps = msa.get_snps(include_ambig=False)

# Consensus sequence
consensus = msa.get_consensus(threshold=0.7, use_ambig_nt=True)

# Reference coordinates
ref_start, ref_end = msa.get_reference_coords()

# Identity matrices
identity_matrix = msa.calc_identity_alignment()      # vs reference or consensus
pairwise_identity = msa.calc_pairwise_identity_matrix(distance_type='ghd')

# Similarity matrices (a variety of substitution matrices available)
similarity = msa.calc_similarity_alignment(matrix_type='BLOSUM65')

# Position weight matrix (DNA/RNA)
pwm = msa.calc_position_matrix(matrix_type='PWM')  # also: 'PFM', 'PPM', 'IC'

# ORFs (DNA/RNA only)
orfs = msa.get_conserved_orfs(min_length=100)
non_overlapping = msa.get_non_overlapping_conserved_orfs(min_length=100)

# Transition/transversion (DNA/RNA)
ts_tv = msa.calc_transition_transversion_score()

# Recovery and frequencies
recovery = msa.calc_percent_recovery()
char_freqs = msa.calc_character_frequencies()

# Alignment statistics
stats = msa.calc_length_stats()  # mean, std, min, max lengths

# Other alignments
reverse_complement = msa.calc_reverse_complement_alignment()  # DNA/RNA
numerical = msa.calc_numerical_alignment()
```

### Plotting

```python
fig, axes = plt.subplots(3, 1, figsize=(14, 10))

# Statistical plot (entropy, gc, or coverage)
draw.stat_plot(msa, ax=axes[0], stat_type='entropy', rolling_average=5)

# Identity alignment
draw.identity_alignment(
    msa, axes[1],
    show_gaps=False,
    show_mismatches=True,
    color_scheme='purine_pyrimidine'
)

# Variant plot (SNPs)
draw.variant_plot(msa, axes[2])
plt.tight_layout()
plt.show()
```

### Available Plot Functions

- `draw.alignment()` - Full sequence visualization
- `draw.identity_alignment()` - Colored by identity to reference
- `draw.similarity_alignment()` - Similarity heatmap
- `draw.stat_plot()` - Entropy, GC, or coverage (+ moving average)
- `draw.variant_plot()` - SNP lolliplots
- `draw.orf_plot()` - ORF annotations (DNA/RNA)
- `draw.annotation_plot()` - Custom feature tracks
- `draw.sequence_logo()` - Sequence conservation
- `draw.consensus_plot()` - Plot consensus sequence

### Data Export

```python
from msaexplorer import explore, export

msa = explore.MSA('alignment.fasta')

# Export SNPs
snps = msa.get_snps()
export.snps(snps, format_type='vcf', path='output.vcf')
export.snps(snps, format_type='tabular', path='output.tsv')

# Export data
export.fasta(msa.alignment, path='output.fasta')
export.stats(msa.calc_entropy(), path='entropy.tsv')
export.character_freq(msa.calc_character_frequencies(), path='freqs.tsv')
export.percent_recovery(msa.calc_percent_recovery(), path='recovery.tsv')

# ORF export (DNA/RNA)
orfs = msa.get_conserved_orfs()
export.orf(orfs, chrom='seq1', path='orfs.gff')
```

### Interactive Web App

Run locally:
```bash
msaexplorer --run
```

Or use online: https://jonas-fuchs.github.io/MSAexplorer/app/

## Supported Formats

**Alignments**: FASTA, CLUSTAL, PHYLIP, STOCKHOLM, NEXUS (auto-detected)

**Annotations**: GenBank (.gb/.gbk), GFF3, BED

**Bio objects**: Bio.Align.MultipleSeqAlignment, Bio.SeqIO.GenBankIterator

## Important Notes

- Reference sequence must exist in alignment
- Annotation sequences must match alignment IDs
- DNA/RNA-only methods raise TypeError on amino acid alignments
- Use `msa.zoom` to focus on regions for faster calculations
- See docstrings: `help(explore.MSA)`, `help(draw.identity_alignment)`, etc.

## License
GNU v.3
See LICENSE file for details.
"""

from importlib.metadata import version, PackageNotFoundError

try:
    __version__ = version("msaexplorer")
except PackageNotFoundError:
    __version__ = "unknown"

