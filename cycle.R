# This script calls the cell cycle.

library(HDF5Array)
data <- readRDS("objects/qc_counts.rds")
hd5mat <- HDF5Array("objects/qc_counts.h5", "neurons") # Need to load both, as HDF5Matrix doesn't have rownames.

library(scran)
mm.pairs <- readRDS(system.file("exdata", "mouse_cycle_markers.rds", package="scran"))
assigned <- cyclone(hd5mat, mm.pairs, gene.names=rownames(data), BPPARAM=MulticoreParam(3))
saveRDS(assigned, file="objects/cycle_output.rds")
