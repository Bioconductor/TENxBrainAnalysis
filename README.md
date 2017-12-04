# Scripts for analyzing the 10X 1.3 million brain cell data

The scripts should be executed in the following order:

- `intro.Rmd`: An introduction, duh.
- `preprocess.Rmd`: Downloading the data and quality control
- `cycle.Rmd`: Cell cycle phase assignment
- `normalize.Rmd`: Calculation of cell-specific size factors
- `variance.Rmd`: Identification of highly variable genes
- `dimred.Rmd`: Dimensionality reduction with randomized PCA 
