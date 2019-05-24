#!/usr/bin/env Rscript

## OUTRGGR ## PART 1: Segmentation algorithm
# Stefano Iantorno ## soap@berkeley.edu
# V.0.1 - May 2019
# Note: verify that all directories are valid prior to running script
# Refer to outrgg.R for more information on each function

source("./outrggr.R")

# Change format as needed between XHMM, CLAMMS, CANOES, and OTHER
# OTHER requires at least 6 columns: "CHR", "START", "END", "INTERVAL", "SAMPLE", "CNV"
# Do quality control and filtering prior to running the algorithm
call.set.1 <- ReadCNVs(read.table("./cnvfile1.xcnv", header = TRUE), format = "XHMM")
call.set.2 <- ReadCNVs(read.table("./cnvfile2.cnv.bed"), format = "CLAMMS")
cnv.list <- list(call.set.1, call.set.2)
names(cnv.list) <- c("XHMM", "CLAMMS")

# This line below is not needed if in OTHER or OUTRGGR format, which should have 8 columns:
# "CHR", "START", "END", "INTERVAL", "SAMPLE", "CNV", "NUM_TARG", "Q_SOME"
cnv.list.ref <- ReformatCNVs(cnv.list, format=names(cnv.list))
cnv.list.ref[[1]] <- bedr.sort.region(cnv.list.ref[[1]], check.chr=F)
cnv.list.ref[[2]] <- bedr.sort.region(cnv.list.ref[[2]], check.chr=F)

# Check if valid regions
test.cond1 <- ifelse(all(is.valid.region(cnv.list.ref[[1]], check.chr=F)), "All regions valid", "Some regions are invalid")
print(test.cond1)
test.cond2 <- ifelse(all(is.valid.region(cnv.list.ref[[2]], check.chr=F)), "All regions valid", "Some regions are invalid")
print(test.cond2)

cnv.list.split <- SplitCNVs(cnv.list.ref, format=names(cnv.list.ref))
overlaps <- Find2Overlaps(cnv.list.split[[1]], cnv.list.split[[2]])

# Set intersect flag to T if interested in intersection rather than union (default)
segments <- CombineSamples(overlaps, merge=F, intersect=T)
write.table(segments, file="./candidate_regions.bed", row.names=F, col.names=F, sep="\t", quote=F)

