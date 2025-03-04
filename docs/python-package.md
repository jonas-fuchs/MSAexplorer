## Installation

MSAexplorer can be used as a simple python extension.

```bash
git clone https://github.com/jonas-fuchs/MSAexplorer
cd MSAexplorer
pip install .
```



### Alignment exploration

#### The MSA class

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

#### The annoation class
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

### Alignment drawing

Simple matplotlib extension to plot alignments. 


**Example:**

```python
import matplotlib.pyplot as plt
from msaexplorer import explore
from msaexplorer import draw

aln = explore.MSA("example_alignments/BoDV.aln", reference_id=None, zoom_range=None)
aln.reference_id = list(aln.alignment.keys())[0]
fig, ax = plt.subplots(nrows=9, height_ratios=[0.2,0.2,0.2,0.2,2,0.2,2,0.2,0.5], sharex=False)

draw.stat_plot(aln, ax[0], "gc", rolling_average=20, line_color="black")
draw.stat_plot(aln, ax[1], stat_type="entropy", rolling_average=1, line_color="indigo")
draw.stat_plot(aln, ax[2], "coverage", rolling_average=1)
draw.stat_plot(aln, ax[3], stat_type="identity", rolling_average=1, line_color="grey")
draw.identity_alignment(aln, ax[4], show_gaps=False, show_mask=True, show_mismatches=True, reference_color='lightsteelblue', show_seq_names=False, show_ambiguities=True, fancy_gaps=True, show_x_label=False, show_legend=True, bbox_to_anchor=(1,1.05))
draw.stat_plot(aln, ax[5], stat_type="similarity", rolling_average=1, line_color="darkblue")
draw.similarity_alignment(aln, ax[6], fancy_gaps=True, show_gaps=True, matrix_type='TRANS', show_cbar=True, cbar_fraction=0.02,  show_x_label=False)
draw.orf_plot(aln, ax[7], cmap='hsv', non_overlapping_orfs=False, show_cbar=True, cbar_fraction=0.2, min_length=150)
draw.annotation_plot(aln, 'example_alignments/AB258389.gff3',  ax[8], feature_to_plot='gene', show_x_label=False)
draw.variant_plot(aln, ax[9], show_x_label=True, show_legend=True, bbox_to_anchor=(1,1.35))

fig.set_size_inches(14, 29)
fig.tight_layout()
plt.show()
```

**Will result in:**
![example](../example_alignments/BoDV.png)
