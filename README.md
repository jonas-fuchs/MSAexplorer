WORK IN PROGRESS... 

## MSAexplorer

_"Explore multiple sequence alignments with a simple python package."_ 

![MSAexplorer](msa_explorer.png)

## Requirements

_Tried to make the requirements as minimal as possible._

- python >= 3.12
- matplotlib >= 3.9
- numpy ~ 2.1.3

## Installation

```bash
git clone https://github.com/jonas-fuchs/MSAexplorer
cd MSAexplorer
pip install .
```

## Documentation

### The MSA class

Allows to read in multiple sequence alignment and compute several statistics.

```python
from msaexplorer import explore

# load MSA
msa = explore.MSA('path-to-msa')

# you can set the zoom range and the reference id
msa.zoom = (0, 1500)  # set zoom range
msa.reference_id = 'your_ref_id'  # set reference id, this is used by different functions (otherwise consensus is used)

# access MSA attributes
msa.alignment  # alignment as dictionary
msa.length  # alignment length
msa.aln_type  # type of the alignment

# explore MSAs, outputs are largely dictionaries, some of the function are only for DNA/RNA
msa.get_snps()  # function to call snps compared to ref (if set) or consensus
msa.get_consensus()  # calculate a consensus sequence
msa.get_conserved_orfs()  # get orfs (including internal) with conserved start and stop position. identity is also calculated.
msa.get_reference_coords()  # get alignment coordinates of the reference
msa.get_non_overlapping_conserved_orfs()  # get orfs that do not overlap, additional a min identity can be set

msa.calc_gc()  # calc gc content
msa.calc_entropy()  # calc shannons normalized entropy
msa.calc_coverage()  # calc coverage
msa.calc_length_stats()  # calculate length statitstics
msa.calc_percent_recovery()  # calculated recovery over reference/consensus sequence
msa.calc_identity_alignment()  # convert the alignment into a identity alignment 
msa.calc_similarity_alignment()  # convert the alignment into a similartity alignment based on different similarity matrices
msa.calc_character_frequencies()  # calculate observed character frequencies
msa.calc_pairwise_identity_matrix()  # calculate identity matrix
msa.calc_reverse_complement_alignment()  # convert alignment sequences to reverese complement
```

### The annoation class
Read in *.gb, *.gff and bed files. All genomic locations are automatically adapted.

```python
from msaexplorer import explore

# read in annotations
aln = explore.MSA('path-to-msa')
annotation = explore.Annotation(aln, 'path-to-annotation')

# access attributes
annotation.ann_type  # type of annotation
annotation.locus  # the locus
annotation.features  # all parsed features matching the corresponsing aln positions

```

## Drawing an alignment with minimal python syntax

```python
import matplotlib.pyplot as plt
from msaexplorer import explore
from msaexplorer import draw

# load alignment
aln = explore.MSA("example_alignments/BoDV.aln", reference_id=None, zoom_range=None)
# define subplots
fig, ax = plt.subplots(nrows=5, height_ratios=[0.5,0.5,0.5,2,1], sharex=False)
# draw different alignment vis
draw.stat_plot(aln, ax[0], "gc", rolling_average=50, line_color="black")
draw.stat_plot(aln, ax[1], stat_type="entropy", rolling_average=50, line_color="indigo")
draw.stat_plot(aln, ax[2], "coverage", rolling_average=1)
draw.identity_alignment(aln, ax[3], show_gaps=True, show_mask=True, show_mismatches=True, reference_color='lightsteelblue', show_seq_names=False, show_ambiguities=True, fancy_gaps=True, show_x_label=False, show_legend=True)
draw.similarity_alignment(aln, ax[4], show_x_label=True)
# format figure
fig.set_size_inches(14, 15)
fig.tight_layout()
plt.show()
```

![example](example_alignments/BoDV.png)