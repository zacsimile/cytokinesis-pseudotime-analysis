# expansion-analysis

Analyze expansion microscopy images of cytokinesis.

## Setup

Clone this repository. Make sure `conda` is installed, ideally through 
[Miniforge](https://github.com/conda-forge/miniforge?tab=readme-ov-file).
Navigate to the folder containing this repository. Then run, in this folder,

```
conda create -n expansion python=3.11
conda activate expansion
pip install -r requirements.txt
```

## Usage

Launch a Jupyter Lab instance in VSCode, another IDE, or through the Miniforge prompt:

```
jupyter lab
```

In Jupyter Lab, open and run the following notebooks in order.

### Notebooks

***step1_find_pseudotime_bins.ipynb*** - Compute image metrics and run PCA on the 
resulting feature space. Divide this space into pseudotime bins. Append the results to 
the original workbooks.

***step2_create_pseudotime_image.ipynb*** - Reconstruct a pseudotime representation
of the midbody, extracted from the original images. The result is an OME-TIFF 
displaying a mean projection, a z-stack, or a radial projection of the midbody over
pseudotime. Color channels represent different target proteins.