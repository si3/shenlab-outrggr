#!/usr/bin/env Rscript
#args = commandArgs(trailingOnly = TRUE) # this will store any arguments to Rscript as entries in string vector args

## outrgg.R script ##
## February 2019 - Columbia University
## Stefano Iantorno (si3 @ github)
## Call using Rscript outrgg.R xcnv1 xcnv2 xcnv3

## Load plyr before dplyr to avoid conflicts
library(plyr)
library(dplyr)

## Make sure that bedtools and bedops are installed
## Check that bedops is set to your number of cores for parallel processing
## Load path to their directories into R session environment
## If not already in your system path

#old_path <- Sys.getenv("PATH")
#Sys.setenv(PATH=paste(old_path, "./.../bin", sep = ":"))
library(bedr)

## Set path to other files
#bamfile_path="./.../bamfile"

## Or read in from args

#if (length(args)==0) {
#  stop("At least two arguments must be supplied.n", call. = FALSE)
#} else if (length(args)==1) {
#  args[2] = "outrggR.batch"
#}

# this loops over all samples
# for big datasets, convert to list of samples
# then run parLapply in parallel package
consensus.list <- vector("list", length(sample.names))
for (SAMPLE in sample.names) {
  i <- which(sample.names == SAMPLE)
  consensus.list[[i]] <- FindConsensus(SAMPLE, counts, cov)
}

# Function to import CNV calls
ReadCNVs <- function(xcnv, format=NULL){
  format<-as.factor(format)
  if (is.null(xcnv)) {stop("Read in CNVs with read.table(), with or without header")}
  if (is.null(format)) {stop("No format given")}
  if (!format %in% c("CANOES", "CLAMMS", "XHMM")) {stop("Please select format between CANOES, CLAMMS, XHMM")}
  # add a function to read.table here
  if (format == "CANOES") {
    col.test <- ifelse(ncol(xcnv)==10, "Correct number of columns for CANOES format", "Incorrect number of columns for CANOES format")
    print(col.test)
    colnames(xcnv) <- c("SAMPLE", "CNV", "INTERVAL", "KB", "CHR", "MID_BP", "TARGETS",
                        "NUM_TARG", "MLCN", "Q_SOME")
  }
  if (format == "CLAMMS") {
    col.test <- ifelse(ncol(xcnv)==18, "Correct number of columns for CLAMMS format", "Incorrect number of columns for CLAMMS format")
    print(col.test)
    colnames(xcnv) <- c("CHR", "START", "END", "INTERVAL", "SAMPLE", "CNV", "MLCN", "NUM_TARG", "Q_SOME", 
                        "Q_EXACT", "Q_LEFT_EXTEND", "LEFT_EXTEND_COORD", "Q_RIGHT_EXTEND", "RIGHT_EXTEND_COORD",
                        "Q_LEFT_CONTRACT", "LEFT_CONTRACT_COORD", "Q_RIGHT_CONTRACT", "RIGHT_CONTRACT_COORD")
  }
  if (format == "XHMM") {
    col.test <- ifelse(ncol(xcnv)==15, "Correct number of columns for XHMM format", "Incorrect number of columns for XHMM format")
    print(col.test)
    colnames(xcnv) <- c("SAMPLE", "CNV", "INTERVAL", "KB", "CHR", "MID_BP", "TARGETS", "NUM_TARG", "Q_EXACT", "Q_SOME", "Q_NON_DIPLOID", "Q_START", "Q_STOP", "MEAN_RD", "MEAN_ORIG_RD")
  }
  return(xcnv)
}

# Function to import count and GC data, nobuild flag doesn't check for formatting if T
ReadData <- function(read.counts, gc, sample.names=NULL, target=NULL, nobuild=F) {
  if (is.null(read.counts)) {stop("Read in count data with read.table(), with or without header")}
  if (is.null(gc)) {stop("Read in GC data as a numeric vector")}
  if (is.null(sample.names)){
    sample.names <- paste("S", seq(4:ncol(read.counts)), sep="")
  }
  if (is.null(target)){
    target <- seq(1, nrow(read.counts))
  }
  if (nobuild){
    return(counts=read.counts)
  }
  names(read.counts) <- c("chromosome", "start", "end", sample.names)
  counts <- cbind(target, gc, read.counts)
  if (ncol(counts) > ncol(read.counts) + 2) {stop("Incorrect number of columns")}
  return(counts)
}

## PART 1 FUNCTIONS ##
# Essential columns include CHR, START, END, INTERVAL, SAMPLE, CNV
# Make sure CHR is character and START, END are numeric

# Function to reformat to BED, accepts CNV list or object (list is slower)
ReformatCNVs <- function(xcnv, format) {
  format <- as.factor(format)
  if (is.null(format)) {stop("Select format between CANOES, CLAMMS, XHMM")}
  if (class(xcnv) == "list"){
    xcnv.list <- as.list(xcnv)
    if (length(xcnv.list) != length(format)) {stop("Insufficient number of formats specified")}
    for (i in seq(1, length(xcnv.list))) {
      xcnv.list[[i]][,3] <- as.character(xcnv.list[[i]][,3])
      # breaks down intervals in CANOES and XHMM, reorders and selects columns
      if (i %in% which(format == "CANOES" | format == "XHMM")) {
        intervals <- t(as.data.frame(strsplit(xcnv.list[[i]][,3], "-")))
        intervals[,1] <- t(as.data.frame(strsplit(intervals[,1], ":")))[,2]
        rownames(intervals) <- NULL
        colnames(intervals) <- c("START", "END")
        xcnv.list[[i]] <- cbind(intervals, xcnv.list[[i]])
      }
      xcnv.list[[i]] <- select(xcnv.list[[i]], "CHR", "START", "END", "INTERVAL", "SAMPLE", "CNV", "NUM_TARG", "Q_SOME")
      xcnv.list[[i]][,2:3] <- apply(xcnv.list[[i]][,2:3], 2, as.integer)
      xcnv.list[[i]][,1] <- as.character(xcnv.list[[i]][,1])
    }
    names(xcnv.list) <- as.character(format)
    return(xcnv.list)
  }
  else {
    if (format == "CANOES" | format == "XHMM") {
      intervals <- t(as.data.frame(strsplit(as.character(xcnv[,3]), "-")))
      intervals[,1] <- t(as.data.frame(strsplit(intervals[,1], ":")))[,2]
      rownames(intervals) <- NULL
      colnames(intervals) <- c("START", "END")
      xcnv <- cbind(intervals, xcnv)
    }
    xcnv <- select(xcnv, "CHR", "START", "END", "INTERVAL", "SAMPLE", "CNV", "NUM_TARG", "Q_SOME")
    xcnv[,2:3] <- apply(xcnv[,2:3], 2, as.integer)
    xcnv[,1] <- as.character(xcnv[,1])
    return(xcnv)
  }
}

# Function to split reformatted xcnv list into sublists of samples
SplitCNVs <- function(xcnv.list, format) {
  if (class(xcnv.list) != "list") {stop("Needs a list of at least 2 different CNV sets")}
  # add function to check for basic colnames
  cnvs.per.sample.list <- vector("list", length(xcnv.list))
  for (i in seq(1, length(xcnv.list))) {
    cnvs.per.sample.list[[i]] <- split(xcnv.list[[i]], xcnv.list[[i]]$SAMPLE)
  }
  names(cnvs.per.sample.list) <- names(xcnv.list)
  # remove samples that have 0 CNVs
  for (i in seq(1, length(cnvs.per.sample.list))) {
    cnvs.per.sample.list[[i]] <- cnvs.per.sample.list[[i]][sapply(cnvs.per.sample.list[[i]], nrow) > 0]
  }
  return(cnvs.per.sample.list)
}

# Accessory functions to get info about overlapping CNV segments
CountOverlaps <- function(overlap.df, n=c("2", "3")) {return(sum(overlap.df$n.overlaps == n))}
GetOverlaps <- function(overlap.df, n=c("2", "3")) {return(filter(overlap.df, overlap.df$n.overlap == n))}

# Function to find overlapping segments from sample sublists, strict flag will output only overlaps (no overhangs) as df
# otherwise will return list of overlaps + overhangs by sample
# make sure input is sorted
Find2Overlaps <- function(xcnv.list, strict = F) {
  if (length(xcnv.list) != 2) {stop("Needs a named list of 2 different CNV sets, split by sample")}
  if (is.null(names(xcnv.list[[1]])) | is.null(names(xcnv.list[[2]]))) {stop("Is the list split by sample? Then specify sample names")}
  x <- which(names(xcnv.list[[1]]) %in% names(xcnv.list[[2]]))
  y <- which(names(xcnv.list[[2]]) %in% names(xcnv.list[[1]]))
  #if (all(names(xcnv.list[[1]][x]) != names(xcnv.list[[2]][y]))) {stop("Sample names do not match")}
  if (length(x) != length(y)) {stop("Sample names do not match")}
  test <- vector("list", length(x))
  test <- mapply(function(a,b) {bedr.join.multiple.region(list(a, b), check.chr=FALSE)}, a=xcnv.list[[1]][x], b=xcnv.list[[2]][y], SIMPLIFY=F)
  # Check if strict flag is set, if it is, only return overlapping segments
  if (strict) {
    test.out <- lapply(test, function(x) GetOverlaps(as.data.frame(x), 2))
    return(test.out)
  }
  return(test) # returns a list of overlapping segments + overhangs organized by sample
}

Find3Overlaps <- function(xcnv.list, strict = F) {
  if (length(xcnv.list) != 3) {stop("Needs a named list of 3 different CNV sets, split by sample")}
  shared <- Reduce(intersect, list(names(xcnv.list[[1]]), names(xcnv.list[[2]]), names(xcnv.list[[3]])))
  x <- which(names(xcnv.list[[1]]) %in% shared)
  y <- which(names(xcnv.list[[2]]) %in% shared)
  z <- which(names(xcnv.list[[3]]) %in% shared)
  if (length(shared) == 0) {stop("Sample names do not match")}
  test <- vector("list", length(shared))
  test <- mapply(function(a,b,c) {bedr.join.multiple.region(list(a,b,c), check.chr=FALSE)}, a=xcnv.list[[1]][x], b=xcnv.list[[2]][y], c=xcnv.list[[3]][z], SIMPLIFY=F)
  # Check if strict flag is set
  if (strict) {
    test.out <- lapply(test, function(x) GetOverlaps(as.data.frame(x), 3))
    return(test.out)
  }
  return(test) # returns list of overlaps + overhangs organized by sample
}

# Function to assign DUP/DEL and Q_SOME score to overlapping segments, accepts list of dataframes
# Only tested using strict flag, make sure sample names have no "." in them
# returns data in user-ready table format
Extract2Overlaps <- function(overlaps, xcnv.list) {
  if (!is.list(xcnv.list)) {
    stop("List of CNVs not split by sample needed as input")
  }
  if (!is.data.frame(xcnv.list[[1]]) | !is.data.frame(xcnv.list[[2]])) {
    stop("List of CNVs in data frame format needed as input")
  }
  if (length(xcnv.list) != 2) {
    stop("List of length 2 needed as input")
  }
  test <- do.call(rbind, overlaps)
  samples <- rownames(test)
  samples <- vapply(strsplit(samples, ".", fixed=T), '[', 1, FUN.VALUE = character(1))
  test <- cbind(test[,1:3], samples)
  colnames(test) <- c("CHR", "START", "END", "SAMPLE")
  overlaps.sorted <- bedr.sort.region(test, check.chr=F)
  xcnv1 <- bedr.sort.region(xcnv.list[[1]], check.chr=F)
  xcnv2 <- bedr.sort.region(xcnv.list[[2]], check.chr=F)
  joined1 <- bedr.join.region(overlaps.sorted, xcnv1, check.chr=F)
  joined2 <- bedr.join.region(overlaps.sorted, xcnv2, check.chr=F)
  # Only consider CNVs called in the same sample
  joined1 <- joined1[joined1[,4]==joined1[,9],]
  joined2 <- joined2[joined2[,4]==joined2[,9],]
  test.out <- cbind(overlaps.sorted, joined1[,c(10,12)], joined2[,c(10,12)])
  colnames(test.out) <- c("CHR", "START", "END", "SAMPLE", "CNV_A", "Q_SOME_A", "CNV_B", "Q_SOME_B")
  return(test.out)
}

Extract3Overlaps <- function(overlaps, xcnv.list) {
  if (!is.list(xcnv.list)) {
    stop("List of CNVs not split by sample needed as input")
  }
  if (!is.data.frame(xcnv.list[[1]]) | !is.data.frame(xcnv.list[[2]]) | !is.data.frame(cnv.list[[3]])) {
    stop("List of CNVs in data frame format needed as input")
  }
  if (length(xcnv.list) != 3) {
    stop("List of length 3 needed as input")
  }
  test <- do.call(rbind, overlaps)
  samples <- rownames(test)
  samples <- vapply(strsplit(samples, ".", fixed=T), '[', 1, FUN.VALUE = character(1))
  test <- cbind(test[,1:3], samples)
  colnames(test) <- c("CHR", "START", "END", "SAMPLE")
  overlaps.sorted <- bedr.sort.region(test, check.chr=F)
  xcnv1 <- bedr.sort.region(xcnv.list[[1]], check.chr=F)
  xcnv2 <- bedr.sort.region(xcnv.list[[2]], check.chr=F)
  xcnv3 <- bedr.sort.region(xcnv.list[[3]], check.chr=F)
  joined1 <- bedr.join.region(overlaps.sorted, xcnv1, check.chr=F)
  joined2 <- bedr.join.region(overlaps.sorted, xcnv2, check.chr=F)
  joined3 <- bedr.join.region(overlaps.sorted, xcnv3, check.chr=F)
  # Only consider CNVs called in the same sample
  joined1 <- joined1[joined1[,4]==joined1[,9],]
  joined2 <- joined2[joined2[,4]==joined2[,9],]
  joined3 <- joined3[joined3[,4]==joined3[,9],]
  test.out <- cbind(overlaps.sorted, joined1[,c(10,12)], joined2[,c(10,12)], joined3[,c(10,12)])
  colnames(test.out) <- c("CHR", "START", "END", "SAMPLE", "CNV_A", "Q_SOME_A", "CNV_B", "Q_SOME_B", "CNV_C", "Q_SOME_C")
  return(test.out)
}

# Function to remove discordant calls, needs list of dataframes, output from Extract2Overlaps
# will return a list of dataframes without discordant calls
RemoveDisc2 <- function(extract.df, xcnv.list) {
  if (length(xcnv.list) != 2) {stop("List of length 2 needed as input")}
  to.remove <- extract.df[extract.df$CNV_A != extract.df$CNV_B,]
  a <- bedr.sort.region(xcnv.list[[1]], check.chr=F)
  b <- bedr.sort.region(xcnv.list[[2]], check.chr=F)
  test.out <- list(a[-which(in.region(a, to.remove, check.chr=F)),], b[-which(in.region(b, to.remove, check.chr=F)),])
  names(test.out) <- names(xcnv.list)
  return(test.out)
}

RemoveDisc3 <- function(extract.df, xcnv.list) {
  if (length(xcnv.list) != 3) {stop("List of length 3 needed as input")}
  to.remove <- extract.df[extract.df$CNV_A != extract.df$CNV_B | extract.df$CNV_A != extract.df$CNV_C | extract.df$CNV_B != extract.df$CNV_C,]
  a <- bedr.sort.region(xcnv.list[[1]], check.chr=F)
  b <- bedr.sort.region(xcnv.list[[2]], check.chr=F)
  c <- bedr.sort.region(xcnv.list[[3]], check.chr=F)
  test.out <- list(a[-which(in.region(a, to.remove, check.chr=F)),], b[-which(in.region(b, to.remove, check.chr=F)),], c[-which(in.region(c, to.remove, check.chr=F)),])
  names(test.out) <- names(xcnv.list)
  return(test.out)
}

# Accessory function to consolidate calls across samples, fetching overhangs
GetClusters <- function(x) {
  test <- cluster.region(x, check.chr=F)
  keep <- names(which(table(test$regionIndex) > 1))
  test.out <- filter(test, regionIndex %in% keep | test$V4 > 1)
  return(test.out)
}

# Function to process output of Find2Overlaps, merge flag will merge intervals across samples, 
# removed strict flag, if by.cluster = F it will consider all intervals
CombineSamples <- function(xcnv.by.sample, merge = F, by.cluster = T) {
  if (!is.list(xcnv.by.sample)) {stop("List of CNVs split by sample needed")}
  if (by.cluster) {
    xcnv.by.sample <- lapply(xcnv.by.sample, GetClusters)
  }
  test <- do.call(rbind, xcnv.by.sample) # faster than ldply, make.row.names=F unused?
  test.sorted <- bedr.sort.region(test[,1:3], check.chr=F) # don't need sample names
  # note that bedr.merge.region has a stratify.by option for merging only in certain groups 
  if ( merge ) {
    test.out <- bedr.merge.region(test.sorted, check.chr=F, list.names=F)
    return(test.out)
  }
  test.out <- bedr(engine="bedops", input=list(test.sorted), method = "partition", params = "", check.chr = F)
  keep <- which(is.valid.region(test.out, check.chr=F))
  test.out <- bedr.sort.region(test.out[keep,], check.chr=F)
  return(test.out)
}


# Functions to call read depth in putative segments
# GetReadDepth # bedtools multicov -bams path/to/bam -bed path/to/bed -q MinMQ
# GetGC # bedtools nuc -fi path/to/fasta -bed path/to/bed (column number 2 after bed entry is %GC)

## PART 2 FUNCTIONS ##

# Calculates correlation matrix for all samples

PrepRefMatrix <- function(counts) {
  if (length(setdiff(names(counts)[1:5], c("target", "chromosome", "start", "end", "gc"))) > 0){
    stop("First five columns of counts matrix must be target, chromosome, start, end, gc")
  }
  counts <- arrange(counts, chromosome, start)
  sample.names <- colnames(counts)[-seq(1,5)]
  mean.counts <- mean(apply(counts[, sample.names], 2, mean))
  counts[, sample.names] <- apply(counts[, sample.names], 2, 
                                  function(x, mean.counts) 
                                    round(x * mean.counts / mean(x)), mean.counts)
  cov <- cor(counts[, sample.names], counts[, sample.names])
  return(cov)
}

# Accessory functions to recalibrate likelihoods in consensus calls, these are called by ModelCopyState

BuildTruthTable <- function(x, states) {
  states.grid <- replicate(nrow(x), states, F)
  states.truth <- expand.grid(states.grid)
  states.truth.table <- apply(states.truth,2,as.character)
  return(states.truth.table)
}

ReplaceProbs <- function(x, states, probs) {
  for (i in 1:ncol(x)) {
    del.i <- which(x[,i] == "DEL")
    dup.i <- which(x[,i] == "DUP")
    dip.i <- which(x[,i] == "DIP")
    x[del.i, i] <- probs[i,"delprob"]
    x[dup.i, i] <- probs[i,"dupprob"]
    x[dip.i, i] <- probs[i,"normalprob"]
  }
  x <- apply(x, 2, as.numeric)
  return(x)
}

# Here below, y is lookup table with sumlogp's, x is the names in merged output
CheckNames <- function(x, y) {
  targets <- unlist(strsplit(x, ","))
  if (length(targets) > 1) {
    if (length(unique(y[targets,"sumlogp1"])) != 1) { 
      #x <- as.matrix(x)
      LR1 <- sum(unique(y[targets, "sumlogpnull"])) - sum(unique(y[targets, "sumlogp1"]))
      LR2 <- sum(unique(y[targets, "sumlogp2"])) - sum(unique(y[targets, "sumlogp1"]))
      x <- c(x, LR1, LR2)
      #sum(unique(y[targets, "sumlogpnull"])) / sum(unique(y[targets, "sumlogp1"]))
      return(x)
    }
    #x <- as.matrix(x)
    LR1 <- unique(y[targets, "sumlogpnull"]) - unique(y[targets, "sumlogp1"])
    LR2 <- unique(y[targets, "sumlogp2"]) - unique(y[targets, "sumlogp1"])
    x <- c(x, LR1, LR2)
    return(x)
  } else {
    #x <- as.matrix(x)
    LR1 <- y[targets, "sumlogpnull"] - y[targets,"sumlogp1"]
    LR2 <- y[targets, "sumlogp2"] - y[targets, "sumlogp1"]
    x <- c(x, LR1, LR2)
    return(x)
  }
}

# Functions to generate likelihoods and output results for each cluster
# ModelCopyState needs segs.by.clust combined with 
ModelCopyState <- function(cluster.index, states, segs.by.clust) {
  # Check that maximum number of intervals is 10 at most
  subcluster <- unique(segs.by.clust[cluster.index, 5])
  probs <- segs.by.clust[cluster.index, 7:9]
  #if (length(subcluster) > 10) stop("Cluster has more 10 segments, but should have no more than 10")
  # Check that you don't have a singleton
  if (length(cluster.index) == 1) {
    copy.state <- apply(segs.by.clust[cluster.index,7:9], 1, which.max)
    state.probs.singleton <- segs.by.clust[cluster.index,7:9]
    sumlogpnull <- state.probs.singleton$normalprob
    sumlogp2 <- apply(state.probs.singleton, 1, function(x) sort(x, partial=length(x)-1)[length(x)-1])
    sumlogp1 <- state.probs.singleton[cbind(seq_along(copy.state), copy.state)]
    copy.state <- names(c(DEL=1, DIP=2, DIP=3)[copy.state])
    results.table <- cbind(segs.by.clust[cluster.index,], copy.state, sumlogp1, sumlogp2, sumlogpnull)
    results.table <- select(results.table, target, copy.state, sumlogp1, sumlogp2, sumlogpnull)
    #results.table$copy.state <- names(c(DEL=1, DIP=2, DIP=3)[copy.state])
    # use select to make matching number of columns then rbind
    return(results.table)
  }
  state.truth <- by(segs.by.clust[cluster.index,], factor(segs.by.clust[cluster.index,5]), BuildTruthTable, states)
  results.list <- vector("list", length(state.truth))
  # Next replace each state with its probability then extract useful variables
  for (i in 1:length(state.truth)) {
    if (ncol(state.truth[[i]]) > 10) stop("Cluster has >10 segments, but should have no more than 10")
    state.probs <- ReplaceProbs(state.truth[[i]], states, filter(segs.by.clust[cluster.index,], subclust==i)[,7:9])
    sumlogp1 <- apply(state.probs, 1, sum)
    sumlogp2 <- sort(sumlogp1, partial=length(sumlogp1)-1)[length(sumlogp1)-1]
    copy.state <- state.truth[[i]][which.max(sumlogp1), ]
    target <- filter(segs.by.clust[cluster.index,], subclust == i)[,6]
    # identify null model and extract null probability
    state.mod <- as.data.frame(cbind(state.truth[[i]], sumlogp1))
    sumlogpnull <- as.numeric(as.character(filter_all(as_tibble(state.mod), all_vars(. != "DEL" & . != "DUP"))$sumlogp1))
    # create results table
    results.table <- cbind.data.frame(target, copy.state, sumlogp1=max(sumlogp1), sumlogp2, sumlogpnull)
    results.list[[i]] <- as.data.frame(results.table)
  }
  results.table <- do.call(rbind, results.list)
  results.table <- as.data.frame(results.table)
  return(results.table)
}

RecalibrateClusters <- function(counts, states=c("DEL","DIP","DUP"), segs.by.clust) {
  # Move all of this below outside of loop
  #segs <- data.frame(chr=as.character(counts[,3]), start=as.numeric(counts[,4]), end=as.numeric(counts[,5]))
  #segs$chr <- as.character(segs$chr)
  #segs$start <- as.numeric(segs$start)
  #segs$end <- as.numeric(segs$end)
  #segs.sort <- bedr.sort.region(segs, check.chr=F)
  #segs.by.clust <- cluster.region(segs.sort, check.chr=F)
  #colnames(segs.by.clust) <- c(colnames(segs), "clust")
  #states <- c("DEL", "DIP", "DUP")
  # Split clusters longer than 10 segments into chunks
  #subcluster.all <- vector("numeric", nrow(segs.by.clust))
  #for (i in unique(segs.by.clust[,4])) {
  #  cluster.i <- which(segs.by.clust[,4] == i)
  #  subcluster <- ceiling(seq_along(segs.by.clust[cluster.i, 4])/10) # can make this 5 for better speed
  #  subcluster.all[cluster.i] <- subcluster
  #}
  #segs.by.clust <- cbind(segs.by.clust, subclust=subcluster.all, state.probs)
  results.list <- vector("list", length(unique(segs.by.clust[,"clust"])))
  # move all of this above outside of loop
  # Run a loop over all cluster
  for (i in 1:length(unique(segs.by.clust[,"clust"]))) {
    cluster.i <- which(segs.by.clust[, "clust"] == i)
    # Apply same code to all clusters,
    #if (length(cluster.i) > 10) {
    #  subcluster.i.list <- split(cluster.i, ceiling(seq_along(cluster.i)/10))
    #  subcluster.states.list <- lapply(subcluster.i.list, ModelCopyState, states = c("DEL", "DIP", "DUP"), segs.by.clust)
    #  # code below needs to wrapped into function and then run with lapply
    #  cluster.states <- do.call(rbind, subcluster.states.list)
    #  #subcluster.results <- lapply(subcluster.i.list, BuildResultsTable, cluster.states, segs.by.clust)
    #  results.table <- do.call(rbind, subcluster.states.list)
    #  results.table <- as.data.frame(results.table)
    #  colnames(results.table) <- c("target", "CNV", "sumlogp1", "sumlogp2", "sumlogpnull")
    #  results.table$target <- as.numeric(as.character(results.table$target))
    #  results.table$CNV <- as.factor(results.table$CNV)
    #  results.table$sumlogp1 <- as.numeric(as.character(results.table$sumlogp1))
    #  results.table$sumlogp2 <- as.numeric(as.character(results.table$sumlogp2))
    #  results.table$sumlogpnull <- as.numeric(as.character(results.table$sumlogpnull))
    #  results.list[[i]] <- results.table
      #print("Cluster split into subclusters")
    #} else if (length(cluster.i) <= 10) { 
      cluster.states <- ModelCopyState(cluster.i, states, segs.by.clust)
      results.table <- as.data.frame(cluster.states)
      colnames(results.table) <- c("target", "CNV", "sumlogp1", "sumlogp2", "sumlogpnull")
      results.table$target <- as.numeric(as.character(results.table$target))
      results.table$CNV <- as.factor(results.table$CNV)
      results.table$sumlogp1 <- as.numeric(as.character(results.table$sumlogp1))
      results.table$sumlogp2 <- as.numeric(as.character(results.table$sumlogp2))
      results.table$sumlogpnull <- as.numeric(as.character(results.table$sumlogpnull))
      results.list[[i]] <- results.table
    #}
  }
  results.table <- do.call(rbind, results.list)
  rownames(results.table) <- NULL
  results.table <- filter(results.table, CNV!="DIP") # consider only non-diploid
  # build lookup table to calculate LR
  lookup <- segs.by.clust[,c(1:3)]
  rownames(lookup) <- segs.by.clust$target
  output.df <- cbind(lookup[as.character(results.table$target),], results.table)
  output.sort.df <- bedr.sort.region(output.df, check.chr=F)
  # Merge neighbouring segments if same CNV type
  merge.output.df <- bedr.merge.region(output.sort.df, stratify.by="CNV", check.chr=F)
  merge.output.df$size <- merge.output.df$end - merge.output.df$start
  merge.output.df <- merge.output.df[merge.output.df$size > 50,] # filter out CNVs smaller than 50bp
  # add code to keep CNV type as column from rownames
  types <- rownames(merge.output.df)
  types <- vapply(strsplit(types, ".", fixed=T), '[', 1, FUN.VALUE = character(1))
  merge.output.df <- data.frame(merge.output.df, type=types)
  # create a lookup table to find the LRs for each target
  lookup <- output.df[,6:8]
  rownames(lookup) <- output.df$target
  merge.out <- sapply(merge.output.df$names, CheckNames, lookup)
  merge.out <- as.data.frame(cbind(merge.out[1,], merge.out[2,], merge.out[3,]))
  rownames(merge.out) <- NULL
  colnames(merge.out) <- c("names", "LR1", "LR2")
  merge.out$names <- as.character(merge.out$names)
  #merge.out <- apply(merge.out, 2, as.character)
  test.out <- left_join(merge.output.df, merge.out)
  test.out$LR1 <- as.numeric(as.character(test.out$LR1))
  test.out$LR2 <- as.numeric(as.character(test.out$LR2))
  return(test.out)
  #return(results.list)
}

# To add above
# script to add SAMPLE column
# script to 

#output <- do.call(rbind, consensus.list)

for (SAMPLE in sample.names) {
  i <- which(sample.names == SAMPLE)
  test.out[[i]] <- data.frame(consensus.list[[i]], sample=SAMPLE)
}

output <- do.call(rbind, consensus.list)
lapply(consensus.list, function(x) {data.frame(x, sample=names(x))})

# Wrapper function to generate consensus calls for all segments, remove.out flag will remove influential points with Cook's D
# Incorporates functions to find clusters and feeds them into ModelCopyState

FindConsensus <- function(sample.name, counts, cov, numrefs=30, maxrows=36000 , remove.out=F) {
  if (!sample.name %in% names(counts)){stop("No column for sample ", sample.name, " in counts matrix")}
  if (length(setdiff(names(counts)[1:5], c("target", "chromosome", "start", "end", "gc"))) > 0){
    stop("First five columns of counts matrix must be target, chromosome, start, end, gc")
  }
  if (numrefs <= 0){
    stop("parameter numrefs must be positive")
  }
  counts <- arrange(counts, chromosome, start)
  sample.names <- colnames(counts)[-seq(1,5)]
  # find mean coverage of probes (mean of means)
  mean.counts <- mean(apply(counts[, sample.names], 2, mean))
  # normalize counts; round so we can use negative binomial
  counts[, sample.names] <- apply(counts[, sample.names], 2, 
                                  function(x, mean.counts) 
                                    round(x * mean.counts / mean(x)), mean.counts)
  reference.samples <- setdiff(sample.names, sample.name)
  covariances <- cov[sample.name, reference.samples]
  reference.samples <- names(sort(covariances, 
                                  decreasing=T)[1:min(numrefs, length(covariances))])
  sample.mean.counts <- mean(counts[, sample.name])
  counts[, reference.samples] <- apply(counts[, reference.samples], 2, 
                                       function(x, sample.mean.counts) 
                                         round(x * sample.mean.counts / 
                                                 mean(x)), sample.mean.counts)
  # Find sample weights
  b <- counts[, sample.name]
  A <- as.matrix(counts[, reference.samples])
  library(nnls)
  all <- nnls(A, b)$x
  est <- matrix(0, nrow=50, ncol=length(reference.samples))
  set.seed(1)
  for (i in 1:50){
    d <- sample(nrow(A), min(500, nrow(A)))
    est[i, ] <- nnls(A[d, ], b[d])$x
  }
  est.weights <- colMeans(est)
  sample.weights <- est.weights / sum(est.weights)
  library(Hmisc)
  # Estimate weighted mean
  counts$mean <- apply(counts[, reference.samples], 
                       1, wtd.mean, sample.weights)
  # Keep non-zero rows
  nonzero.rows <- counts$mean > 0
  nonzero.rows.df <- data.frame(target=counts$target, 
                                nonzero.rows=nonzero.rows)
  counts <- counts[nonzero.rows, ]
  # Estimate weighted variance
  counts$var <- apply(counts[, reference.samples], 1, wtd.var, sample.weights, normwt=T) # takes a minute
  set.seed(1)
  counts.subset <- counts[sample(nrow(counts), min(maxrows, nrow(counts))), ]
  library(mgcv)
  counts.subset$var[counts.subset$var==0] <- 0.1 
  # Remove outliers if remove.out is T by fitting a glm approximation and calculating Cook's D
  #if (remove.out) {
  #  fit.test <- glm(var ~ mean + gc, family="Gamma", data=counts.subset)
  #  dist <- cooks.distance(fit.test)
  #  remove <- which(dist > 4/maxrows)
  #  counts.subset <- counts.subset[-remove,]
  #}
  # fitting a GAM to estimate parameters, remove outliers if remove.out is T
  fit <- gam(var ~ s(mean) + s(gc), family=Gamma(link=log), data=counts.subset)
  if (remove.out) {
    dist <- cooks.distance(fit)
    remove <- which(dist > 4/maxrows)
    counts.subset <- counts.subset[-remove,]
    fit <- gam(var ~ s(mean) + s(gc), family=Gamma(link=log), data=counts.subset)
  }
  # choose between predicted value, weighted variance, or Poisson
  v.estimate <- pmax(predict(fit, counts, type="response"), counts$var, 
                     counts$mean * 1.01)
  var.estimate <- cbind(counts$target, v.estimate)
  test.counts <- counts[,sample.name]
  target.means <- counts$mean 
  var.estimate <- v.estimate
  targets <- counts[, "target"]
  num.targets <- length(counts[,sample.name])
  # calculate the means for the deletion, normal and duplication states
  state.target.means <- t(apply(data.frame(x=counts$mean), 1, function(x) c(x*1/2, x, x*3/2)))
  # calculate the expected size (given the predicted variance)
  size <- (counts$mean ^ 2) / (v.estimate - counts$mean)
  state.probs <- matrix(NA, num.targets, 4)
  colnames(state.probs) <- c("target", "delprob", "normalprob", "dupprob")
  # calculate the probabilities of del, dip, dup at each interval given the read count
  size.del <- size / 2
  size.dup <- size * 3 / 2
  state.probs[, "delprob"] <- dnbinom(counts[, sample.name], mu=state.target.means[, 1], size=size.del, log=T)
  state.probs[, "normalprob"] <- dnbinom(counts[, sample.name], mu=state.target.means[, 2], size=size, log=T)
  state.probs[, "dupprob"] <- dnbinom(counts[, sample.name], mu=state.target.means[, 3], size=size.dup, log=T)
  state.probs[, "target"] <- counts[, "target"]
  # some values may be infinite as a result of extreme read count
  # these must be removed
  row.all.inf <- which(apply(state.probs, 1, function(x){all(is.infinite(x))}))
  if (length(row.all.inf) > 0){
    for (i in row.all.inf){
      if (test.counts[i] >= state.target.means[i, 3]){
        state.probs[i, 2:4] <- c(-Inf, -Inf, -0.01)
      }
      else if (test.counts[i] <= state.target.means[i, 1]){
        state.probs[i, 2:4] <- c(-0.01, -Inf, -Inf)
      }
      else state.probs[i, 2:4] <- c(-Inf, -0.01, -Inf)
    }
  }
  # now generate clusters and truth table for each cluster
  segs <- data.frame(chr=as.character(counts[,3]), start=as.numeric(counts[,4]), end=as.numeric(counts[,5]))
  segs$chr <- as.character(segs$chr)
  segs$start <- as.numeric(segs$start)
  segs$end <- as.numeric(segs$end)
  segs.sort <- bedr.sort.region(segs, check.chr=F)
  segs.by.clust <- cluster.region(segs.sort, check.chr=F)
  colnames(segs.by.clust) <- c(colnames(segs), "clust")
  states <- c("DEL", "DIP", "DUP")
  # Split clusters longer than 10 segments into chunks
  subcluster.all <- vector("numeric", nrow(segs.by.clust))
  for (i in unique(segs.by.clust[,4])) {
    cluster.i <- which(segs.by.clust[,4] == i)
    subcluster <- ceiling(seq_along(segs.by.clust[cluster.i, 4])/10) # can make this 5 for better speed
    subcluster.all[cluster.i] <- subcluster
  }
  segs.by.clust <- cbind(segs.by.clust, subclust=subcluster.all, state.probs)
  test.out <- RecalibrateClusters(counts, states, segs.by.clust)
  return(test.out)
  #return(state.probs)
}

Phred <- function(LR) {
  score <- round(min(99, -10 * log10(exp(LR))))
  return(score)
}

#sapply(test$LR1, Phred)
