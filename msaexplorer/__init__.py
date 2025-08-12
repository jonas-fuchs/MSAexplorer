r"""
# What is MSAexplorer?

MSAexplorer allows the analysis and straight forward plotting of multiple sequence alignments.
Its focus is to act as a simple python3 extension or shiny app with minimal dependencies and syntax. It's easy
to set up and highly customizable.

# Installation

#### Via pip (recommended)
```bash
pip install msaexplorer # or
pip install msaexplorer[process]  # additionally installs pyfamsa and pytrimal (not required, but optional in the app)
```

#### From this repo
```bash
git clone https://github.com/jonas-fuchs/MSAexplorer
cd MSAexplorer
pip install . # or
pip install .[process]
```

# Usage as a shiny application

The current version of the app is also deployed to [GitHub pages](https://jonas-fuchs.github.io/MSAexplorer/app/). This application is serverless, and all
computation runs through your browser. There is no need to install anything. Enjoy the app!

However, you can also run it yourself or host it however you like!

#### Running the app
```bash
msaexplorer --run
```
Now just follow the link provided in your terminal.

#### Exporting as a static site
```bash
pip install shinylive
git clone https://github.com/jonas-fuchs/MSAexplorer
cd MSAexplorer
shinylive export ./ site/  # you should now have a new 'site' folder with the app
```

# Usage as a python3 package

If you only want to use the MSAexplorer package without the shiny app, you can install it as follows:

```bash
pip install msaexplorer
```

## Analysis

The `explore` module lets you load an alignment file and analyze it.

```python
'''
a small example on how to use the MSAexplorer package
'''

from msaexplorer import explore

# load MSA
msa = explore.MSA('example_alignments/DNA.fasta')
annotation = explore.Annotation(msa, 'example_alignments/DNA_RNA.gff3')

# you can set the zoom range and the reference id if needed
msa.zoom = (0, 1500)
msa.reference_id = 'your_ref_id'

# access functions on what to compute on the MSA
msa.calc_pairwise_identity_matrix()
```

Importantly, multiple sequence alignments should have the format:

```
>Seq1
ATCGATCGATCGATCG
>Seq2
ATCGATCGATCGATCG
>Seq3
ATCGATCGATCGATCG
```

Additionally, you can also read in an annotation in `bed`, `gff3` or `gb` format and connect them to the MSA. Importantly,
the sequence identifier has to be part of the alignment. All genomic locations are then automatically adapted to the
alignment.

## Plotting

The plotting `draw` module has several predefined functions to draw alignments.

```python
'''
an example demonstrating how to plot multiple sequence alignments
'''
# import all packages
import matplotlib.pyplot as plt
from msaexplorer import explore
from msaexplorer import draw

#  load alignment
aln = explore.MSA("example_alignments/DNA.fasta", reference_id=None, zoom_range=None)
# set reference to e.g. the first sequence
aln.reference_id = list(aln.alignment.keys())[0]

fig, ax = plt.subplots(nrows=2, height_ratios=[0.2, 2], sharex=False)

draw.stat_plot(
    aln,
    ax[0],
    stat_type="entropy",
    rolling_average=1,
    line_color="indigo"
)

draw.identity_alignment(
    aln,
    ax[1],
    show_gaps=False,
    show_mask=True,
    show_mismatches=True,
    reference_color='lightsteelblue',
    color_scheme='purine_pyrimidine',
    show_seq_names=False,
    show_ambiguities=True,
    fancy_gaps=True,
    show_x_label=False,
    show_legend=True,
    bbox_to_anchor=(1,1.05)
)

plt.show()
```
"""

from importlib.metadata import version, PackageNotFoundError

try:
    __version__ = version("bamdash")
except PackageNotFoundError:
    __version__ = "unknown"