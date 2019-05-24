# OUTRGGR: an integrated probabilistic framework for consensus calling of CNVs from whole-exome sequencing data #

OUTRGGR-CANOES (outlier-adjusted regression via GAM for genotyping in R of CNVs with an arbitrary number of exome samples), or simply OUTRGGR, is a platform to assess read depth support for CNV calls made by multiple CNV detection tools starting from WES data. It provides an integrated approach to consensus CNV calling, boundary resolution, and recalibration of genotype likelihood scores. 

## Software requirements: ##

OUTRGGR is implemented as a series of functions in R script, and relies heavily on the bedr package for set operations. For bedr to work properly, bedtools, bedops, and tabix need to be installed on your system and their directories need to be in your R environment’s path. To manually add the path to your bedtools, bedops, and tabix executables, type the following in your R session or within your R script:

`old_path <- Sys.getenv("PATH")`

`Sys.setenv(PATH=paste(old_path, "path/to/bin", sep = ":"))`

Please refer to bedr package materials for more detailed instructions on how to install it and configure it. Other R package requirements include: plyr, dplyr, nnls, Hmisc, and mgcv. These can be installed by running the following command.

`install.packages(c(“bedr”, “plyr”, “dplyr”, “nnls”, “Hmisc”, “mgcv”))`

If working with a large dataset, we recommend using the parallel package in R.

## Overview: ##

OUTRGGR is designed to run in two separate steps, although they can be easily combined in a single pipeline. In the first step (the **“segmentation”** step), the script accepts CNV calls in BED format as input, finds overlapping segments, and returns a BED file with either only the overlapping segments (option “strict”), all overlapping segments plus any overhanging segments (option “intersect”), or all possible segments (option “union”). The default is to consider all segments, or the union of original CNV calls. To maximize specificity, the intersection plus any overhangs can be used instead. By definition, "union" segments will include "intersect" segments.

The user then estimates read depth in these candidate CNV intervals using the original BAM or CRAM files for the study samples. In the second step (the **“modelling”** step), the interval coordinates, read depth, and GC information for all samples (in CANOES format) is used as input. The algorithm performs a series of normalization, filtering, regression, and parameter estimation steps to build models of all possible copy state changes at each consecutive interval, and outputs the most likely model across each candidate CNV region, along with a likelihood score.

A few built-in functions are provided to read in XHMM, CLAMMS, or CANOES formatted CNV calls. CNV files from any other CNV caller should be converted into OUTRGGR format manually. We require the following columns in each BED file: 

* CHR – chromosome number
* START – coordinate of start position
* END – coordinate of end positon
* INTERVAL – interval coordinates (in chr:start-end format)
* SAMPLE – sample identifier
* CNV – type of CNV (DEL or DUP)

> NOTE: for bedr to recognize this object as a BED file once read into R, the CHR column needs to be of class “character”, and the START and END columns need to be of class “numeric”. 

The current version of the algorithm supports comparisons of either two or three different calling platforms. Given the computational complexity of integrating more than 2 datasets, with unclear benefit in performance, only the combination of 2 sets has been tested. It is important to choose methods that have complementary strengths and weaknesses in terms of their performance, as the resulting CNV consensus calls will have properties that reflect a balance of both input datasets. Input CNV calls should have undergone QC, as the final quality of the output consensus calls will depend on the quality of the input calls.

> NOTE: all the built-in functions that are part of the bedr package should work as described in the package documentation, and can be used to customize the script. 

Additional filtering should be performed on the consensus calls based on the Phred scaled LR1 score, which is a measure of the likelihood of the alternative non diploid model, and the LR2 score, which is a measure of the likelihood of the most likely model over the second most likely model, and thus indicates how well the data supports the overall CNV boundaries.

## Workflow: ##

The segmentation algorithm can be run as a standalone R script as follows:

`nohup Rscript Run_Segmentation.R > nohup.segs.log &`

Please refer to the script for modification of variables and relevant directories based on your system.

Since bedtools does not support CRAM format at this time, the calculation of read depth in candidate segments should be done with either samtools (for CRAM files) or bedtools (for BAM files): 

`samtools bedcov path/to/bed path/to/cram`

`bedtools multicov -bams path/to/bam -bed path/to/bed -q MinMQ`

Since `samtools bedcov` calculates the sum of the read depth at each individual base pair position, and `bedtools multicov` calculates the absolute read count over the whole interval, both are expected to follow the negative binomial distribution, unlike tools that calculate the mean read depth over the interval of interest. We use a minimum mapping quality (`MinMQ`) of 20. This step should be parallelized using GNU `parallel` or `xargs` if working with a large dataset.

To estimate GC content, we use the reference fasta file and run the following:

`bedtools nuc -fi path/to/fasta -bed path/to/bed`

The resulting text file will have %GC content in the second column after each BED entry. After calculating read depth in the candidate segments and estimating GC content, we can then model copy number states over all candidate segments:

`nohup Rscript Run_Modelling.R > nohup.model.log &`

Depending on the sample size of your dataset, you may need to parallelize your script. We use the parallel package in R to speed up computation. Instead of running the previous line, a parallelized script can be used:

`nohup Rscript Run_Modelling_Parallel.R > nohup.model.par.log &`

Processing about 2000 samples on a machine with ~30 cores in this manner took 2 hours (excluding read depth calculation time, which can be done either within R or invoking external bash scripts according to user preference). The algorithm can thus be easily integrated into existing CNV calling pipelines.

## References: ##

* bedr package - https://cran.r-project.org/web/packages/bedr/vignettes/Using-bedr.html
* bedtools - https://bedtools.readthedocs.io/en/latest/
* CANOES - http://www.columbia.edu/~ys2411/canoes/
