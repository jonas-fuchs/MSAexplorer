![Logo](app_src/img/logo.svg)

[![language](https://img.shields.io/badge/python-%3E3.11-green)](https://www.python.org/)
[![License: GPL v3](https://img.shields.io/github/license/jonas-fuchs/bamdash)](https://www.gnu.org/licenses/gpl-3.0)

_**Some personal thoughts**, because let's face it: Generating beautiful alignment plots is always a hassle. There are a lot of alignment viewers and cool tools out there.
For all the Geneious users who now think there ain't a problem... While this software has its perks (well it should, as it has the costs of a minivan), I am waiting since 10 years that you can plot a 
half decent alignment with it. So if you feel the same way, you are probably also switching between tools which is time-consuming and not really satisfying. This is why 
I wanted something that is primarily developed for generating high quality publication ready figures while retaining simplicity. Most of all I wanted 
flexibility. I had the vision to generate nice plots of alignments no matter if I just want to plot parts or the whole alignment. 
Moreover, I wanted something to quickly generate and plot statistics without having to rely on multiple tools. So this is why I developed MSAexplorer and I really hope its
something for you! I think it turned out pretty cool but I let you be the judge of that :smiley:._

### MSAexplorer is available as a python package and also a standalone app to analyse multiple sequence alignments and generate publication ready figures. 
**Want to just use MSAexplorer and generate publication ready figures? 
The curently stable version of the MSAexplorer app is hosted on  [github pages](https://jonas-fuchs.github.io/MSAexplorer/app).**


### Features of MSAexplorer as an app
- :white_check_mark: The app runs solely in your browser. No need to install anything, just might take a few seconds to load.
- :white_check_mark: Use the app offline (after loading it).
- :white_check_mark: Analyse alignments on your smartphone or tablet.
- :white_check_mark: Download alignment statistics (e.g. entropy, SNPs, coverage, consensus, ORFs and more).
- :white_check_mark: Annotate the alignment by additionally reading in gb, gff or bed files.
- :white_check_mark: Flexibility to customize plots and colors.
- :white_check_mark: Easily export the plot as pdf.
- :white_check_mark: Generate plots of the whole alignment as well as just parts of it.
- :white_check_mark: Publication ready figures with just a few clicks.

| ![](readme_assets/upload_tab.png) | ![](readme_assets/plot_tab.png) | ![](readme_assets/plot2_tab.png) | ![](readme_assets/analysis_tab.png) |
|-----------------------------------|---------------------------------|----------------------------------|-------------------------------------|


### Features of MSAexplorer as a python package ([full documentation](https://jonas-fuchs.github.io/MSAexplorer/docs/msaexplorer.html))
- :white_check_mark: Access MSAexplorer as a python package
- :white_check_mark: Maximum flexibility for the plotting and analysis features while retaining minimal syntax.
- :white_check_mark: Integrates seamlessly with matplotlib.
- :white_check_mark: Minimal requirements.

```python
### Minimal analysis example ###

from msaexplorer import explore, export

# load the alignment
aln = explore.MSA('example_alignments/DNA.fasta')
print(aln.aln_type)  # print alignment type
print(aln.length)  # print alignment length

# adjust what you want to look at
aln.reference_id = 'AB032031.1 Borna disease virus 1 genomic RNA, complete genome'  # set a reference if needed
aln.zoom = (0, 1000)  # set a zoom range

# now print for example all snps in that zoom range compared to the reference id
snps = aln.get_snps(include_ambig=True)
print(snps)
# and then save to file
export.snps(snps, path='my_path/snps.vcf', format_type='vcf')

# see documentation for full usage
```


```python
### Two minimal plotting examples ###

from msaexplorer import explore, draw
import matplotlib.pyplot as plt

# Example 1
draw.identity_alignment('example_alignments/DNA.fasta')
plt.show()

# Example 2
aln = explore.MSA('example_alignments/DNA.fasta')
# adjust zoom levels (for example to also plot sequence text)
aln.zoom = (0,60)
plt.figure(figsize=(10,12))  # adjust so the sequence text fits in your figure well
draw.identity_alignment(aln, show_identity_sequence=True)
plt.show()

```

```python
### Extended plotting example  ####

import matplotlib.pyplot as plt
from msaexplorer import explore
from msaexplorer import draw

#  load alignment
aln = explore.MSA("example_alignments/DNA.fasta", reference_id=None, zoom_range=None)
# set reference to first sequence
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

## Requirements

`python >= python 3.11`

The requirements for the python package have been kept minimal:
- `matplotlib>=3.8`
- `numpy>=2.0`

If you want to use the MSAexplorer app then you will additionally need to install:
- `shiny>=1.3`
- `shinywidgets>=0.5.2`
- `plotly>=5.23`
- `pyfamsa>=0.5.3`

## Installation
### MSAexplorer python package
```bash
# via pip
pip install msaexplorer
# from this repo
git clone https://github.com/jonas-fuchs/MSAexplorer
cd MSAexplorer
pip install .
```
### MSAexplorer app
#### Local installation
Run it locally on your machine with:
````bash
git clone https://github.com/jonas-fuchs/MSAexplorer
cd MSAexplorer
pip install msaexplorer  # installs the msaexplorer package
pip install -r requirements.txt  # installs app dependencies
shiny run app.py
````
#### Hosting MSAexplorer
If you want to host MSAexplorer e.g. for your group, you can export the app as a static html with a few easy steps:
```bash
# install shinylive for exporting
pip install shinylive
git clone https://github.com/jonas-fuchs/MSAexplorer
cd MSAexplorer
shinylive export ./ site/  # you should now have a new 'site' folder with the app
```