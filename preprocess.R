# This script converts the 1M Neuron data set into a HDF5Matrix.

path <- "1M_neurons_filtered_gene_bc_matrices_h5.h5"

library(TENxGenomics)
tenx.se <- tenxSummarizedExperiment(path)
tenx <- assay(tenx.se, withDimnames=FALSE)

##############################################
# Calculating margin summaries.

margin.summary <- function(x, nrow) {
    ## > str(x)
    ## List of 3
    ##  $ ridx : num [1:20548381] 8 9 17 39 52 63 118 119 123 182 ...
    ##  $ cidx : num [1:20548381] 1 1 1 1 1 1 1 1 1 1 ...
    ##  $ value: int [1:20548381] 1 1 2 2 1 7 2 1 1 1 ...

    ## rows: summarize all rows, whether in current sample or not.
    ridx <- structure(                  # quick 'factor'
        x$ridx, .Label=as.character(seq_len(nrow)), class="factor"
    )
    rowdf <- data.frame(
        ridx = seq_len(nrow),
        n = tabulate(x$ridx, nrow),
        sum = vapply(split(x$value, ridx), sum, numeric(1), USE.NAMES=FALSE)
    )

    ## columns: summarized cells (complete) in current sample
    ucidx <- unique(x$cidx)
    x$cidx <- match(x$cidx, ucidx)
    coldf <- data.frame(
        cidx = ucidx,
        n = tabulate(x$cidx, length(ucidx)),
        sum = vapply(split(x$value, x$cidx), sum, numeric(1), USE.NAMES=FALSE)
    )

    list(rowdf = rowdf, coldf = coldf)
}

library(BiocParallel)
register(MulticoreParam(progressbar=TRUE))
result <- tenxiterate(tenx, margin.summary, nrow = nrow(tenx))
rows <- Reduce(function(x, y) {
    idx <- c("n", "sum")
    x[, idx] <- x[, idx] + y[, idx]
    x
}, lapply(result, `[[`, 1))
cols <- do.call("rbind", lapply(result, `[[`, 2))
colnames(cols) <- c("Cell", "Ngenes", "Libsize")

# Saving all cell-based data to file.
dir.create("objects", showWarnings=FALSE)
write.table(file="objects/cell_data.txt", cbind(cols, colData(tenx.se)[,-(1:2)]), 
            sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

##############################################
# Doing some basic quality control on the cells.

lib.sizes <- cols$Libsize
nexprs <- cols$Ngenes
batch <- tenx.se$Library

library(scater)
keep <- !isOutlier(lib.sizes, log=TRUE, nmad=3, batch=batch, type="lower") &
        !isOutlier(nexprs, log=TRUE, nmad=3, batch=batch, type="lower") 
saveRDS(keep, file="objects/qc_keep.rds")

# Subsetting the TENxMatrix object.
tenxmat <- TENxMatrix(path)[,keep]

##############################################
# Choosing chunk dimensions.

#ncells <- ncol(tenxmat)
#ngenes <- nrow(tenxmat)
#cache.size <- ncells*sqrt(ngenes)
#chunk.nrow <- ceiling(cache.size/ncells)
#chunk.ncol <- ceiling(cache.size/ngenes)

# Saving to file as column-based chunks for now.
library(HDF5Array)
options(DelayedArray.block.size=2e8)
out <- writeHDF5Array(tenxmat, file="objects/qc_counts.h5", name="neurons", chunk_dim=c(1e4, 1))
rownames(out) <- rownames(tenxmat)
saveRDS(out, file="objects/qc_counts.rds")

