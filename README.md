![Logo](docs/logo.svg)


Introducing MSAexplorer, a simple python package to analyse and plot multiple sequence alignments.
The main goal of the plotting features are to grasp the whole alignment which is particular useful 
for large alignments and where the goal is to generate a publication ready overviews. Moreover, MSAexplorers
simplifies reading and analyses of multiple sequence alignments with minimal python syntax.

### Usage

MSAexplorer can be used in different ways, depending on your needs.

* [**Shiny app:**](docs/shiny-app.md) Access MSAexplorer's plotting features directly through a [shiny UI](https://shiny.posit.co/py/). Plots are highly customizable. [For a live demo click me!](https://foxlabs.shinyapps.io/msaexplorer/)
* [**Command line tool:**](docs/standalone.md) Access MSAexplorer through the terminal. (Coming soon...)
* [**Python extension:**](docs/python-package.md) Directly access the packages [API](https://jonas-fuchs.github.io/MSAexplorer/) to analyse and plot MSA with simple python syntax.


### Requirements
* python >= 3.12

### Quick Installation
```bash
git clone https://github.com/jonas-fuchs/MSAexplorer
cd MSAexplorer
pip install .
```