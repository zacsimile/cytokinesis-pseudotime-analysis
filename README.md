# expansion-analysis

This is code for analyzing expansion images of cytokinesis.

## Setup

Clone this repository. Make sure you have `conda`, ideally install it through 
[Miniforge](https://github.com/conda-forge/miniforge?tab=readme-ov-file).
Naviate to the folder containing this repository. Then run, in this folder,

```
conda create -n expansion python=3.11
conda activate expansion
pip install -r requirements.txt
```

## Usage

Launch a Jupyter Lab in VSCode, another IDE, or through the Miniforge prompt:

```
jupyter lab
```

In Jupyer Lab, open the correct notebook to run.

## Notebook descriptions

***images_over_pseudotime.ipynb*** - Reconstruct a pseudotime representation of a single experiment (e.g. Septin2-GFP). The result will be an OME TIFF displaying mean projections across time points (e.g. defined by RS, CS, etc. or by PCA sorting) and color channels representing the different structures.