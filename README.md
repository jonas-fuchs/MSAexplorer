![Logo](app_src/img/logo.svg)

Introducing MSAexplorer, a simple python package and [shiny application](https://shiny.posit.co/py/) to analyse and plot multiple sequence alignments.
The main goal of the plotting features are to grasp the whole alignment and generate publication ready pdfs. Moreover, MSAexplorer
simplifies reading and analyses of multiple sequence alignments with minimal python syntax.

The serverless shiny app is automatically deployed to gihub pages via [shinylive](https://shiny.posit.co/py/docs/shinylive.html) and includes the current version of the master branch. All computation runs in your browser. MSAexplorer can handle any alignment size. However, the computation time increases with larger alignments. If your alignment is similar to my [examples](example_alignments/DNA.fasta), plotting will be near to instant on any modern cpu.

### Requirements

`python >= python 3.11`

## Quick Installation (python package)

```bash
git clone https://github.com/jonas-fuchs/MSAexplorer
cd MSAexplorer
pip install .
```
## [Shiny app (Click me)](https://jonas-fuchs.github.io/MSAexplorer/shiny)
Runs in your browser - no need to install anything, just might take a few seconds to load.

Want to install it locally?
````bash
git clone https://github.com/jonas-fuchs/MSAexplorer
cd MSAexplorer
pip install .  # installs the msaexplorer package
pip install -r requirements.txt  # installs shiny dependencies
shiny run app.py
````


![png](readme_assets/shiny-app.png)

## [Python package documentation](https://jonas-fuchs.github.io/MSAexplorer/docs/msaexplorer.html)

### [Exploration](https://jonas-fuchs.github.io/MSAexplorer/docs/msaexplorer/explore.html) 

### [Plotting](https://jonas-fuchs.github.io/MSAexplorer/docs/msaexplorer/draw.html) 

#### [Mayhem example _(lets just go nuts )_ :bomb: ](https://jonas-fuchs.github.io/MSAexplorer/docs/msaexplorer.html#plotting)

![png](readme_assets/BoDV.png)
