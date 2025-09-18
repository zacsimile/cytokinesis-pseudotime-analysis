# cytokinesis-pseudotime-analysis

Split 3D expansion microscopy images of cytokinesis into pseudotime windows based on image features. Reconstruct average mean projection, radial projection, and z-stack images using these pseudotime windows.

## Requirements

- Windows 10+, Mac 13.7+ or Linux (Ubuntu 20.04+). 8 GB RAM. *NOTE: Only tested on Windows and Mac.*
- [Miniforge](https://github.com/conda-forge/miniforge?tab=readme-ov-file) or an equivalent `conda` environment manager.

## Setup

Clone this repository. Navigate to the folder containing this repository. Then run, in this folder,

```
conda create -n expansion python=3.11
conda activate expansion
pip install -r requirements.txt
```

## Usage

### targets.yaml

All notebooks and scripts rely on this file to direct them to data locations. For each protein target of interest, add an entry to this file containing

* The "key", `ALIXm` in this example. This is the most commonly used reference to this target in the raw file names in the `image_directory` (see below).
* `image_directory`: the location of the raw images
* `workbook`: path to an XLSX file containing information about the the images
* `workbook_header_row`: the row in the workbook where the column names are loated
* `workbook_sheet_name`: which sheet contains data about the raw images
* `alias` *(optional)*: if the raw file names refer to this target by more than one name, enter alternatives here. All aliases will be mapped to the "key" in the results for consistent naming of this target.

```
ALIXm:
  image_directory: /path/to/ALIX m
  workbook: /path/to/20250114_ALIXm.xlsx
  workbook_header_row: 0
  workbook_sheet_name: Tabelle1
  alias:
  - ALIX
  - ALIXM
```

### Notebooks

Launch a Jupyter Lab instance in VSCode, another IDE, or through the Miniforge prompt:

```
jupyter lab
```

In Jupyter Lab, open and run the following notebooks in order.

***step1_find_pseudotime_bins.ipynb*** - Compute image metrics and run PCA on the 
resulting feature space. Divide this space into pseudotime bins. Append the results to 
the original workbooks.

***step2_create_pseudotime_image.ipynb*** - Using the results of `step1_find_pseudotime_bins.ipynb`, reconstruct a pseudotime representation
of the midbody, extracted from the original images. The result is an OME-TIFF 
displaying a mean projection, a z-stack, or a radial projection of the midbody over
pseudotime. Color channels represent different target proteins.

***looped_images_over_pseudotime copy.ipynb*** - Provided with pre-calculated pseudotime windows, reconstruct a pseudotime representation
of the midbody, extracted from the original images. The result is an OME-TIFF 
displaying a mean projection, a z-stack, or a radial projection of the midbody over
pseudotime. Color channels represent different target proteins.

### Scripts

Run from a command line (Windows) or Terminal (Mac):

```
conda activate expansion
python looped_images_over_pseudotime.py
```

***looped_images_over_pseudotime.py*** - Performs the same function as `looped_images_over_pseudotime copy.ipynb`, but from a Python script.

## Related datasets

A subset of the data that this package was developed for, is stored as a [Zenodo repository](https://doi.org/10.5281/zenodo.17143816). It contains a limited number of example data for the targets actin and MKLP1. It also contains the workbooks corresponding to the data.
To use the notebooks with the provided data subsets, targets.yaml  has to be updated with the targets, name and location of the notebooks and the directory of the downloaded data. 
In `looped_images_over_pseudotime copy.ipynb` the following changes have to be made:

```
time_key = "Stage"

time_order = ["RC", "CS", "RS", "SM", "BA", "A"]

time_do_not_fit = ["RC", "CS", "A"]

desired_channel_order = ["MTs", "septin", "DAPI", "MKLP1", "actin"]
```

## Related repositories 

The code to generate the table for the merged metrics option, which includes features extracted from Cellprofiler, can be found in the [cytokinesis-feature-analysis](https://github.com/AG-Ewers/cytokinesis-feature-analysis) repository.
