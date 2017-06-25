library(HDF5Array)
hd5mat <- HDF5Array("objects/norm_exprs.h5", "neurons")
cell.data <- read.table("objects/cell_data.txt", header=TRUE)
kept <- readRDS("objects/qc_keep.rds")

# Only considering genes that are expressed at a decent level.
library(scater)
ab <- calcAverage(hd5mat, size.factors=cell.data$Libsize[kept])
keep <- ab >= 0.01

# Blocking on the library of origin.
lib.source <- factor(cell.data$Library[kept])
design <- model.matrix(~0 + lib.source)

# Estimating the variance.
library(scran)
fit <- trendVar(hd5mat, trend="semiloess", design=design, subset.row=keep)
dec <- decomposeVar(fit=fit)
rownames(dec) <- which(keep)
write.table(file="objects/hvg_output.txt", dec, sep="\t", quote=FALSE, col.names=NA)
