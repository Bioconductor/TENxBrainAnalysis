# Scripts for analyzing the 10X 1.3 million brain cell data

The scripts should be executed in the following order:

- `intro.Rmd`: An introduction, duh.
- `preprocess.Rmd`: Downloading the data and quality control
- `cycle.Rmd`: Cell cycle phase assignment
- `normalize.Rmd`: Calculation of cell-specific size factors
- `variance.Rmd`: Identification of highly variable genes
- `dimred.Rmd`: Dimensionality reduction with randomized PCA 

Various output objects will be saved to `objects/`.
A few of these objects are currently hosted at https://drive.google.com/open?id=1_0WbmJ2BriLKlyKEf1Bbb8K0_NwD9rw-.
Note that `sce.rds` does not contain the actual counts or normalized expression values, and requires something like this:

```r
library(TENxBrainData)
tenx <- TENxBrainData()
sce <- readRDS("sce.rds")
counts(sce) <- counts(tenx) # overwrite inbuilt absolute path
sce <- normalize(sce) # generate normalized expression values
```

The `pics/make_pics.R` scripts will generate the figures used in the paper.
