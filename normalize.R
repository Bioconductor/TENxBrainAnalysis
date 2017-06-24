# This script performs normalization.

library(HDF5Array)
hd5mat <- HDF5Array("objects/qc_counts.h5", "neurons")
cell.data <- read.table("cell_data.txt", header=TRUE)
kept <- readRDS("objects/qc_keep.rds")

# Selecting reasonably high-abundance genes for normalization.
library(scater)
ab <- calcAverage(hd5mat, size.factors=cell.data$Libsize[kept])
keep <- ab >= 0.1

# Splitting by batch, and also within batch 
# (the QR decomposition forms a full matrix and requires
# too much memory for large numbers of cells).
lib.source <- factor(cell.data$Library[kept])
by.lib <- split(rep(1:5, length.out=length(lib.source)), lib.source)
for (x in names(by.lib)) {
    by.lib[[x]] <- paste0(x, ".", sort(by.lib[[x]]))       
}
mod.source <- unlist(by.lib)

#current <- hd5mat[,mod.source=="1.1"]
#ab <- rowMeans(current) 
#keep <- ab >= 0.1
#current <- current[keep,]
#
#xcurrent <- realize(current, "HDF5Array")
#sf <- scran:::.computeSumFactors(xcurrent)

# Computing size factors.
library(scran)
sfs <- computeSumFactors(hd5mat, cluster=mod.source, subset.row=keep)
saveRDS(assigned, file="objects/norm_output.rds")

# Computing normalized log-expression values.
options(DelayedArray.block.size=2e8)
log.norm <- log2(t(t(hd5mat)/sfs) + 1)

ncells <- ncol(hd5mat)
ngenes <- nrow(hd5mat)
cache.size <- ncells*sqrt(ngenes)
chunk.nrow <- ceiling(cache.size/ncells)
chunk.ncol <- ceiling(cache.size/ngenes)

writeHDF5Array(log.norm, file="objects/norm_exprs.h5", name="neurons", chunk=c(chunk.nrow, chunk.ncol))
