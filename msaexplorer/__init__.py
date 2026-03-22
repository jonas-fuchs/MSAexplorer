r"""
# MSAexplorer

MSAexplorer is a lightweight toolkit to **explore**, **plot**, and **export** multiple sequence alignments.
The API is organized into three modules:

- `explore`: parse alignments/annotations and compute statistics.
- `draw`: create matplotlib plots directly from `MSA` objects or file paths.
- `export`: serialize computed results (SNPs, stats, FASTA, ORFs, recovery tables).

All examples below are tested with files from `example_alignments/`.

## Quick Start (`explore`)

```python
from pathlib import Path
from msaexplorer import explore

base = Path('example_alignments')
aln = explore.MSA(str(base / 'DNA.fasta'))

# set reference and zoom window (start inclusive, end exclusive)
aln.reference_id = aln.sequence_ids[0]
aln.zoom = (0, 300)

print(aln.aln_type)        # 'DNA'
print(len(aln), aln.length)
print(aln.get_reference_coords())
```

### Biopython Interoperability

```python
from pathlib import Path
from Bio import AlignIO, SeqIO
from msaexplorer import explore

base = Path('example_alignments')

# alignment from Bio.Align.MultipleSeqAlignment
bio_aln = AlignIO.read(str(base / 'DNA.fasta'), 'fasta')
aln = explore.MSA(bio_aln, reference_id=bio_aln[0].id, zoom_range=(0, 200))

# annotation from Bio.SeqIO GenBank iterator
gb_iter = SeqIO.parse(str(base / 'DNA_RNA.gb'), 'genbank')
ann = explore.Annotation(aln, gb_iter)
print(ann.ann_type, list(ann.features.keys())[:3])
```

### Working with Downstream Dataclasses

```python
from pathlib import Path
from msaexplorer import explore

aln = explore.MSA(str(Path('example_alignments') / 'DNA.fasta'), zoom_range=(0, 300))
aln.reference_id = aln.sequence_ids[0]

entropy = aln.calc_entropy()                              # AlignmentStats
length_stats = aln.calc_length_stats()                    # LengthStats
dist_to_ref = aln.calc_pairwise_distance_to_reference()   # PairwiseDistance
variants = aln.get_snps()                                 # VariantCollection
orfs = aln.get_non_overlapping_conserved_orfs(90)         # OrfCollection

print(entropy.stat_name, entropy.positions[:3], entropy.values[:3])
print(length_stats.mean_length, length_stats.std_length)
print(dist_to_ref.reference_id, dist_to_ref.sequence_ids[:2], dist_to_ref.distances[:2])
print(variants.chrom, len(variants))
print(orfs.keys()[:3])
```

## Plotting (`draw`)

All plotting functions return a matplotlib `Axes`.
You can either:

1. pass an existing axis (best for multi-panel figures), or
2. pass only an alignment/path for a one-liner plot.

### One-Liner Examples (path input)

```python
from msaexplorer import draw
import matplotlib.pyplot as plt

draw.identity_alignment('example_alignments/DNA.fasta')
plt.show()
draw.stat_plot('example_alignments/DNA.fasta', stat_type='entropy', rolling_average=5)
plt.show()
```

### Multi-Panel Figure with Main Plot Types

```python
# import necessary packages
from pathlib import Path
import matplotlib.pyplot as plt
from msaexplorer import explore, draw

# Example 1
# load DNA alignment
base = Path('example_alignments')
aln = explore.MSA(str(base / 'DNA.fasta'), zoom_range=(0, 100))

# ini figure
fig, ax = plt.subplots(nrows=3, ncols=1, figsize=(14, 15), height_ratios=[1, 1, 10])
# plot
draw.stat_plot(aln, ax=ax[0], stat_type='coverage', rolling_average=1, show_title=True)
draw.sequence_logo(aln, ax=ax[1], plot_type='logo')
draw.identity_alignment(aln, ax=ax[2], show_consensus=True, show_seq_names=True, fancy_gaps=True, show_legend=True, color_scheme='standard', show_identity_sequence=True)

plt.tight_layout()
plt.show()

# Example 2
# load AA alignment
aln = explore.MSA(str(base / 'AS.fasta'))
# ini figure
fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(15, 5), height_ratios=[1, 2])

draw.stat_plot(aln, ax=ax[0], stat_type='entropy', rolling_average=5, show_title=True)
draw.identity_alignment(aln, ax=ax[1], show_consensus=True, show_seq_names=True, fancy_gaps=True, show_legend=True, color_scheme='hydrophobicity')

plt.tight_layout()
plt.show()
```

## Export (`export`)

```python
from pathlib import Path
from msaexplorer import explore, export

aln = explore.MSA(str(Path('example_alignments') / 'DNA.fasta'), zoom_range=(0, 300))
aln.reference_id = aln.sequence_ids[0]

variants = aln.get_snps()
entropy = aln.calc_entropy()
consensus = aln.get_consensus()

# return as strings
vcf_text = export.snps(variants, format_type='vcf')
stats_text = export.stats(entropy)
fasta_text = export.fasta(consensus, header='consensus_dna')

# or write to disk (path argument)
# export.snps(variants, format_type='tabular', path='results/snps')
# export.stats(entropy, path='results/entropy.tsv')
```

## Notes

- `zoom` can be reset with `aln.zoom = None`.
- If no `reference_id` is set, methods that need a reference use the consensus.
- For API-level details, inspect module/class docstrings in `explore`, `draw`, and `export`.

"""

from importlib.metadata import version, PackageNotFoundError

try:
    __version__ = version("msaexplorer")
except PackageNotFoundError:
    __version__ = "unknown"
