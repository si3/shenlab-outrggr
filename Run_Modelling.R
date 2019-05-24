#!/usr/bin/env Rscript

## OUTRGGR ## PART 2: Modelling algorithm
# Stefano Iantorno ## soap@berkeley.edu
# V.0.1 - May 2019
# Note: verify that all directories are valid prior to running script
# Refer to outrgg.R for more information on each function
# If large sample size, use parallelized script

source("./outrggr.R")

# Requires GC file, read depth file, and sample name list
gc.df <- read.table("./candidate_regions_nuc.bed") # output from bedtools nuc
rd.df <- read.table("./counts.RD.txt") # output from samtools or bedtools, in CANOES format
samples <- read.table("./samples.txt") # list of sample names
counts <- ReadData(rd.df, gc.df$V5, sample.names=as.character(samples[,1]))

cov <- PrepRefMatrix(counts)
sample.names <- colnames(counts)[6:ncol(counts)]
consensus.list <- vector("list", length(sample.names))

# Set remove.out flag to T to improve GAM fit
for (SAMPLE in sample.names) {
  i <- which(sample.names == SAMPLE)
  consensus.list[[i]] <- FindConsensus(SAMPLE, counts, cov, remove.out=F)
}

names(consensus.list) <- sample.names
results.out <- do.call(rbind, consensus.list)
samples <- rownames(results.out)
samples <- vapply(strsplit(samples, ".", fixed=T), '[', 1, FUN.VALUE = character(1))
results.out <- cbind(results.out, sample=samples)
write.table(results.out, "./cnvs.bed", row.names=F, col.names=F, sep="\t", quote=F)

