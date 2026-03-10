r"""
# What is MSAexplorer?

MSAexplorer is a comprehensive Python package for analyzing and visualizing multiple sequence alignments (MSAs).
It combines powerful statistical analysis, publication-quality plotting, and flexible data export in a simple,
dependency-light package. Perfect for both interactive use and integration into bioinformatics pipelines.

## Key Features

- **Multiple input formats**: FASTA, CLUSTAL, PHYLIP, STOCKHOLM, NEXUS, or direct Biopython objects
- **Annotation support**: GenBank, GFF3, BED formats with automatic coordinate mapping
- **Rich statistics**: Entropy, GC content, coverage, pairwise identity, ORF detection, SNP analysis, and more
- **Publication-ready plots**: Identity matrices, statistical plots, similarity heatmaps, and custom visualizations
- **Flexible export**: Generate reports, export alignments in multiple formats, process with external tools
- **Interactive web app**: Shiny-based interface with no installation required
- **Biopython integration**: Direct support for Bio.SeqIO and Bio.AlignIO objects for seamless workflows

## Installation

### Via pip (recommended)
```bash
pip install msaexplorer
# or with optional sequence processing tools:
pip install msaexplorer[process]  # adds pyfamsa and pytrimal support
```

### From source
```bash
git clone https://github.com/jonas-fuchs/MSAexplorer
cd MSAexplorer
pip install .  # or pip install .[process]
```

## Quick Start

### As a Web Application

Launch the interactive shiny app:
```bash
msaexplorer --run
```

Or use the web version (no installation): [GitHub Pages](https://jonas-fuchs.github.io/MSAexplorer/app/)

Export as a static site:
```bash
pip install shinylive
shinylive export ./ site/
```

### As a Python Package

#### Basic Analysis
```python
from msaexplorer import explore

# Load alignment (supports FASTA, CLUSTAL, PHYLIP, STOCKHOLM, NEXUS)
msa = explore.MSA('alignment.fasta')

# Or from a Biopython object
from Bio import AlignIO
bio_alignment = AlignIO.read('alignment.fasta', 'fasta')
msa = explore.MSA(bio_alignment)

# Get basic statistics
print(f"Sequences: {len(msa.alignment)}")
print(f"Length: {msa.length} bp")
print(f"Type: {msa.aln_type}")  # DNA, RNA, or AA

# Compute statistics
entropy = msa.calc_entropy()
gc_content = msa.calc_gc()
coverage = msa.calc_coverage()
pairwise_identity = msa.calc_pairwise_identity_matrix()
```

#### With Annotations
```python
from msaexplorer import explore

# Load alignment
msa = explore.MSA('alignment.fasta')

# Load annotation (GenBank, GFF3, or BED format)
annotation = explore.Annotation(msa, 'annotation.gb')

# Or from a Biopython GenBank iterator
from Bio import SeqIO
gb_iterator = SeqIO.parse('annotation.gb', 'genbank')
annotation = explore.Annotation(msa, gb_iterator)

# Access annotation features
print(annotation.features.keys())
print(f"Annotation type: {annotation.ann_type}")
```

#### Advanced Analysis
```python
from msaexplorer import explore

msa = explore.MSA('alignment.fasta')

# Set reference and zoom range
msa.reference_id = 'seq1'
msa.zoom = (0, 1000)

# Statistical analyses
length_stats = msa.calc_length_stats()
snps = msa.get_snps(include_ambig=False)
consensus = msa.get_consensus(threshold=0.7, use_ambig_nt=True)

# For nucleotide alignments
identity_matrix = msa.calc_identity_alignment()
similarity_matrix = msa.calc_similarity_alignment()
position_matrix = msa.calc_position_matrix(matrix_type='PWM')

# For DNA/RNA alignments
reverse_complement = msa.calc_reverse_complement_alignment()
conserved_orfs = msa.get_conserved_orfs(min_length=100)
ts_tv_score = msa.calc_transition_transversion_score()

# Recovery statistics
recovery = msa.calc_percent_recovery()
char_frequencies = msa.calc_character_frequencies()
```

### Plotting and Visualization

#### Basic Identity Plot
```python
import matplotlib.pyplot as plt
from msaexplorer import explore, draw

# Load and prepare alignment
aln = explore.MSA('alignment.fasta', reference_id='seq1')

# Create identity visualization
fig, ax = plt.subplots(figsize=(14, 8))
draw.identity_alignment(
    aln,
    ax,
    show_gaps=False,
    show_mask=True,
    show_mismatches=True,
    color_scheme='purine_pyrimidine',
    show_seq_names=True,
    show_legend=True
)
plt.tight_layout()
plt.show()
```

#### Statistical Plots
```python
import matplotlib.pyplot as plt
from msaexplorer import explore, draw

aln = explore.MSA('alignment.fasta')

fig, axes = plt.subplots(nrows=4, figsize=(14, 10), sharex=True)

# Entropy plot
draw.stat_plot(aln, axes[0], stat_type='entropy', rolling_average=5)
axes[0].set_ylabel('Entropy')

# GC content plot
draw.stat_plot(aln, axes[1], stat_type='gc', rolling_average=5)
axes[1].set_ylabel('GC Content')

# Coverage plot
draw.stat_plot(aln, axes[2], stat_type='coverage', rolling_average=5)
axes[2].set_ylabel('Coverage')

# SNP frequency
snps = aln.get_snps()
draw.snp_plot(aln, axes[3], snps=snps)
axes[3].set_ylabel('SNP Count')

plt.tight_layout()
plt.show()
```

#### Similarity Heatmaps
```python
import matplotlib.pyplot as plt
from msaexplorer import explore, draw

aln = explore.MSA('alignment.fasta', reference_id='seq1')

fig, ax = plt.subplots(figsize=(10, 8))
draw.similarity_alignment(
    aln,
    ax,
    matrix_type='similarity'  # or 'identity', 'numerical'
)
plt.tight_layout()
plt.show()
```

#### Comparison Plots
```python
import matplotlib.pyplot as plt
from msaexplorer import explore, draw

aln1 = explore.MSA('alignment1.fasta', reference_id='seq1')
aln2 = explore.MSA('alignment2.fasta', reference_id='seq1')

fig, (ax1, ax2) = plt.subplots(ncols=2, figsize=(16, 6))

draw.identity_alignment(aln1, ax1)
draw.identity_alignment(aln2, ax2)

ax1.set_title('Alignment 1')
ax2.set_title('Alignment 2')

plt.tight_layout()
plt.show()
```

### Data Export and Processing

```python
from msaexplorer import explore, export

msa = explore.MSA('alignment.fasta')

# Export to different formats
export.to_fasta(msa, 'output.fasta')
export.to_phylip(msa, 'output.phy')
export.to_nexus(msa, 'output.nex')

# Export statistics
stats = msa.calc_length_stats()
# Process statistics as needed

# Get alignment data
alignment_dict = msa.alignment  # dict[seq_id: sequence]
consensus = msa.get_consensus()
identity_matrix = msa.calc_pairwise_identity_matrix()
```

## Input Format Support

### Alignments
Automatically detects format or specify explicitly:
- **FASTA** - Most common, widely supported
- **CLUSTAL** - From CLUSTAL/MUSCLE alignments
- **PHYLIP** - PHYLIP sequential format
- **STOCKHOLM** - Pfam/HMMER format
- **NEXUS** - PAUP/MrBayes format
- **Bio.Align.MultipleSeqAlignment** - Direct Biopython objects

### Annotations
Supports standard bioinformatics formats:
- **GenBank** (.gb, .gbk) - Feature-rich annotation format
- **GFF3** - General Feature Format v3
- **BED** - Browser Extensible Data format
- **Bio.SeqIO.GenBankIterator** - Direct Biopython iterators

## Configuration and Customization

### Set Reference and Zoom
```python
from msaexplorer import explore

msa = explore.MSA('alignment.fasta')

# Set reference sequence for identity calculations
msa.reference_id = 'my_reference'

# Focus on a region
msa.zoom = (100, 500)  # or just start position
msa.zoom = 100  # equals (100, alignment_end)

# Reset zoom
msa.zoom = None
```

### Color Schemes
Available color schemes for plotting:
- `purine_pyrimidine` - Distinguishes chemical properties
- `nucleotide` - Standard ATCG coloring
- `clustalx` - ClustalX color scheme
- `taylor` - Taylor amino acid coloring
- And more...

## Workflow Examples

### Pipeline: Analyze and Visualize
```python
from msaexplorer import explore, draw
import matplotlib.pyplot as plt

# Load data
msa = explore.MSA('alignment.fasta')
annotation = explore.Annotation(msa, 'annotation.gff3')

# Set parameters
msa.reference_id = list(msa.alignment.keys())[0]
msa.zoom = (0, 2000)

# Compute statistics
stats = {
    'entropy': msa.calc_entropy(),
    'gc': msa.calc_gc(),
    'coverage': msa.calc_coverage(),
    'snps': msa.get_snps(),
    'identity': msa.calc_pairwise_identity_matrix()
}

# Create comprehensive figure
fig = plt.figure(figsize=(16, 12))
gs = fig.add_gridspec(3, 2, hspace=0.3, wspace=0.3)

# Plot 1: Entropy
ax1 = fig.add_subplot(gs[0, :])
draw.stat_plot(msa, ax1, stat_type='entropy')
ax1.set_ylabel('Entropy')

# Plot 2: Identity alignment
ax2 = fig.add_subplot(gs[1:, :])
draw.identity_alignment(msa, ax2, show_seq_names=True)

plt.tight_layout()
plt.show()
```

### Pipeline: Comparative Genomics
```python
from msaexplorer import explore

# Load multiple alignments
alignments = {
    'gene_a': explore.MSA('gene_a.fasta'),
    'gene_b': explore.MSA('gene_b.fasta'),
    'gene_c': explore.MSA('gene_c.fasta'),
}

# Compute comparative statistics
results = {}
for gene_name, msa in alignments.items():
    results[gene_name] = {
        'length_stats': msa.calc_length_stats(),
        'entropy': msa.calc_entropy(),
        'snps': msa.get_snps(),
        'pairwise_identity': msa.calc_pairwise_identity_matrix()
    }

# Use results for downstream analysis
for gene_name, stats in results.items():
    print(f"{gene_name}:")
    print(f"  Mean length: {stats['length_stats']['mean length']:.0f} bp")
    print(f"  SNP count: {len(stats['snps']['POS'])}")
```

## Documentation

For detailed API documentation, use Python's built-in help:
```python
from msaexplorer import explore, draw, export

help(explore.MSA)
help(explore.Annotation)
help(draw.identity_alignment)
help(export.to_fasta)
```

Or view the full documentation at the [GitHub repository](https://github.com/jonas-fuchs/MSAexplorer).

## Citation

If you use MSAexplorer in your research, please cite:
```
[Citation information to be added]
```

## License

MIT License - See LICENSE file for details.
"""

from importlib.metadata import version, PackageNotFoundError

try:
    __version__ = version("msaexplorer")
except PackageNotFoundError:
    __version__ = "unknown"
