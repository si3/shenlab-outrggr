#!/usr/bin/env Rscript

## OUTRGGR ## PART 2: Modelling algorithm
# Stefano Iantorno ## soap@berkeley.edu
# V.0.1 - May 2019
# Note: verify that all directories are valid prior to running script
# Refer to outrgg.R for more information on each function
# Using parallel package - if errors, run the non-parallel version for debugging

source("./outrggr.R")
library(parallel)
no_cores <- detectCores() - 1

# Requires GC file, read depth file, and sample name list
gc.df <- read.table("./candidate_regions_nuc.bed") # output from bedtools nuc
rd.df <- read.table("./counts.RD.txt") # output from samtools or bedtools, in CANOES format
samples <- read.table("./samples.txt") # list of sample names
counts <- ReadData(rd.df, gc.df$V5, sample.names=as.character(samples[,1]))

cov <- PrepRefMatrix(counts)
sample.names <- colnames(counts)[6:ncol(counts)]
consensus.list <- vector("list", length(sample.names))
cl <- makeCluster(no_cores, type="FORK")

# Set remove.out flag to T to improve GAM fit
consensus.list <- parLapply(cl, sample.names, function(x) FindConsensus(x, counts, cov, remove.out=F))
stopCluster(cl)

names(consensus.list) <- sample.names
results.out <- do.call(rbind, consensus.list)
samples <- rownames(results.out)
samples <- vapply(strsplit(samples, ".", fixed=T), '[', 1, FUN.VALUE = character(1))
results.out <- cbind(results.out, sample=samples)
write.table(results.out, "./cnvs.bed", row.names=F, col.names=F, sep="\t", quote=F)

