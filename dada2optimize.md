dada2-optimize
================
Michael G. LaMontagne
1/7/2022

## R Markdown

This R Markdown document has scripts that can be used to optimize the
selection of parameters in DADA2. The goal is to select parameters that
maximize the number of ASVs shared between technical replicates and
maximize the confidence of matching reads against reference sequences.

This script requires pairs of fastq files (for and rev) generated from
four samples. It helps these four samples consist of two duplicate
samples. In this demo we run four samples (G6S, JA25I, JA25II, SH13Z).
G6S and SH13Z are a pair of technical replicates and JA25I and JA25II
are a pair.

## Preprocessing. This step runs in a conda environment. Use bioconda to install fastq-tools and then check the installation. This step is run in Linux systems.

``` bash
#conda create --name fastqTools fastq-tools
#conda activate fastqTools
# random reads from a fastq file
#fastq-sample --help
```

## Subsample 1,000 reads from each pair of fastq files following this example. This step is also run in Linux.

``` bash
#fastq-sample G6S_S36_L001_R1_001.fastq G6S_S36_L001_R2_001.fastq -n 1000 #-s 50 -o G6S_test_
```

## Install packages

``` r
library(dada2); packageVersion("dada2")
```

    ## Loading required package: Rcpp

    ## [1] '1.20.0'

``` r
library(philentropy); packageVersion("philentropy")
```

    ## [1] '0.5.0'

``` r
library(DECIPHER); packageVersion("DECIPHER")
```

    ## Loading required package: Biostrings

    ## Loading required package: BiocGenerics

    ## Loading required package: parallel

    ## 
    ## Attaching package: 'BiocGenerics'

    ## The following objects are masked from 'package:parallel':
    ## 
    ##     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
    ##     clusterExport, clusterMap, parApply, parCapply, parLapply,
    ##     parLapplyLB, parRapply, parSapply, parSapplyLB

    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs

    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, append, as.data.frame, basename, cbind, colnames,
    ##     dirname, do.call, duplicated, eval, evalq, Filter, Find, get, grep,
    ##     grepl, intersect, is.unsorted, lapply, Map, mapply, match, mget,
    ##     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
    ##     rbind, Reduce, rownames, sapply, setdiff, sort, table, tapply,
    ##     union, unique, unsplit, which.max, which.min

    ## Loading required package: S4Vectors

    ## Loading required package: stats4

    ## 
    ## Attaching package: 'S4Vectors'

    ## The following objects are masked from 'package:base':
    ## 
    ##     expand.grid, I, unname

    ## Loading required package: IRanges

    ## 
    ## Attaching package: 'IRanges'

    ## The following object is masked from 'package:philentropy':
    ## 
    ##     distance

    ## The following object is masked from 'package:grDevices':
    ## 
    ##     windows

    ## Loading required package: XVector

    ## Loading required package: GenomeInfoDb

    ## 
    ## Attaching package: 'Biostrings'

    ## The following object is masked from 'package:base':
    ## 
    ##     strsplit

    ## Loading required package: RSQLite

    ## [1] '2.20.0'

``` r
sessionInfo()
```

    ## R version 4.1.1 (2021-08-10)
    ## Platform: x86_64-w64-mingw32/x64 (64-bit)
    ## Running under: Windows 10 x64 (build 22000)
    ## 
    ## Matrix products: default
    ## 
    ## locale:
    ## [1] LC_COLLATE=English_United States.1252 
    ## [2] LC_CTYPE=English_United States.1252   
    ## [3] LC_MONETARY=English_United States.1252
    ## [4] LC_NUMERIC=C                          
    ## [5] LC_TIME=English_United States.1252    
    ## 
    ## attached base packages:
    ## [1] stats4    parallel  stats     graphics  grDevices utils     datasets 
    ## [8] methods   base     
    ## 
    ## other attached packages:
    ##  [1] DECIPHER_2.20.0     RSQLite_2.2.8       Biostrings_2.60.2  
    ##  [4] GenomeInfoDb_1.28.2 XVector_0.32.0      IRanges_2.26.0     
    ##  [7] S4Vectors_0.30.0    BiocGenerics_0.38.0 philentropy_0.5.0  
    ## [10] dada2_1.20.0        Rcpp_1.0.7         
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] lattice_0.20-44             png_0.1-7                  
    ##  [3] Rsamtools_2.8.0             assertthat_0.2.1           
    ##  [5] digest_0.6.27               utf8_1.2.2                 
    ##  [7] plyr_1.8.6                  R6_2.5.1                   
    ##  [9] ShortRead_1.50.0            evaluate_0.14              
    ## [11] ggplot2_3.3.5               pillar_1.6.4               
    ## [13] zlibbioc_1.38.0             rlang_0.4.11               
    ## [15] blob_1.2.2                  Matrix_1.3-4               
    ## [17] rmarkdown_2.11              BiocParallel_1.26.2        
    ## [19] stringr_1.4.0               bit_4.0.4                  
    ## [21] RCurl_1.98-1.5              munsell_0.5.0              
    ## [23] DelayedArray_0.18.0         compiler_4.1.1             
    ## [25] xfun_0.25                   pkgconfig_2.0.3            
    ## [27] htmltools_0.5.2             tidyselect_1.1.1           
    ## [29] SummarizedExperiment_1.22.0 tibble_3.1.4               
    ## [31] GenomeInfoDbData_1.2.6      matrixStats_0.61.0         
    ## [33] fansi_0.5.0                 crayon_1.4.1               
    ## [35] dplyr_1.0.7                 GenomicAlignments_1.28.0   
    ## [37] bitops_1.0-7                grid_4.1.1                 
    ## [39] gtable_0.3.0                lifecycle_1.0.1            
    ## [41] DBI_1.1.1                   magrittr_2.0.1             
    ## [43] scales_1.1.1                RcppParallel_5.1.4         
    ## [45] cachem_1.0.6                stringi_1.7.5              
    ## [47] reshape2_1.4.4              hwriter_1.3.2              
    ## [49] latticeExtra_0.6-29         ellipsis_0.3.2             
    ## [51] generics_0.1.1              vctrs_0.3.8                
    ## [53] RColorBrewer_1.1-2          tools_4.1.1                
    ## [55] bit64_4.0.5                 Biobase_2.52.0             
    ## [57] glue_1.4.2                  purrr_0.3.4                
    ## [59] jpeg_0.1-9                  MatrixGenerics_1.4.3       
    ## [61] fastmap_1.1.0               yaml_2.2.1                 
    ## [63] colorspace_2.0-2            GenomicRanges_1.44.0       
    ## [65] memoise_2.0.0               knitr_1.36

## Create blank variables

``` r
# 
x <- c(0.1:0.5)
D <- matrix(x, nrow=18, ncol=1)
trackRun <- matrix(x, nrow=6, ncol=1)
rownames(D) <- c("try", "t1", "t2", "t3", "Rich", "G12", "H12", "GH", "NoChim", "M12", "BootClass", "BootGenus", "input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim" )
# reset counters
try <- c(0)
M12 <- c(0)
D
```

    ##           [,1]
    ## try        0.1
    ## t1         0.1
    ## t2         0.1
    ## t3         0.1
    ## Rich       0.1
    ## G12        0.1
    ## H12        0.1
    ## GH         0.1
    ## NoChim     0.1
    ## M12        0.1
    ## BootClass  0.1
    ## BootGenus  0.1
    ## input      0.1
    ## filtered   0.1
    ## denoisedF  0.1
    ## denoisedR  0.1
    ## merged     0.1
    ## nonchim    0.1

## Create path to folder with four fastq pairs

``` r
path <- "C:/Users/mglam/Dropbox/BIOINFO/Rscripts/HHdada/UNTreads/UNTtest" 
# CHANGE ME to the directory containing the fastq files after unzipping.
list.files(path)
```

    ## [1] "filtered"                         "G6S_test_L001_R1_001.fastq.gz"   
    ## [3] "G6S_test_L001_R2_001.fastq.gz"    "JA25I_test_L001_R1_001.fastq.gz" 
    ## [5] "JA25I_test_L001_R2_001.fastq.gz"  "JA25II_test_L001_R1_001.fastq.gz"
    ## [7] "JA25II_test_L001_R2_001.fastq.gz" "SH13Z_test_L001_R1_001.fastq.gz" 
    ## [9] "SH13Z_test_L001_R2_001.fastq.gz"

## Forward and reverse fastq filenames must have format: SAMPLENAME\_R1\_001.fastq and SAMPLENAME\_R2\_001.fastq

``` r
fnFs <- sort(list.files(path, pattern="_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
## plot some quality profiles from forward
plotQualityProfile(fnFs[1:2])
```

    ## Warning: `guides(<scale> = FALSE)` is deprecated. Please use `guides(<scale> =
    ## "none")` instead.

![](dada2optimize_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

``` r
## and reverse
plotQualityProfile(fnRs[1:2])
```

    ## Warning: `guides(<scale> = FALSE)` is deprecated. Please use `guides(<scale> =
    ## "none")` instead.

![](dada2optimize_files/figure-gfm/unnamed-chunk-6-2.png)<!-- -->

## Place filtered files in filtered/ subdirectory

``` r
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
head(filtFs)
```

    ## [1] "C:/Users/mglam/Dropbox/BIOINFO/Rscripts/HHdada/UNTreads/UNTtest/filtered/G6S_F_filt.fastq.gz"   
    ## [2] "C:/Users/mglam/Dropbox/BIOINFO/Rscripts/HHdada/UNTreads/UNTtest/filtered/JA25I_F_filt.fastq.gz" 
    ## [3] "C:/Users/mglam/Dropbox/BIOINFO/Rscripts/HHdada/UNTreads/UNTtest/filtered/JA25II_F_filt.fastq.gz"
    ## [4] "C:/Users/mglam/Dropbox/BIOINFO/Rscripts/HHdada/UNTreads/UNTtest/filtered/SH13Z_F_filt.fastq.gz"

## Select range of parameters based on inspection of above quality plots and run optimization loop.

``` r
repeat{
t1 <- (sample(180:240, 1))
t2 <- (sample(140:200, 1))
t3 <- (sample(16:20, 1))
## count trys
try <- try +1
## standard filtering 
# On Windows set multithread=FALSE
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, trimLeft=t3, truncLen=c(t1,t2),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE) 
head(out)
## learn error rates
errR <- learnErrors(filtRs, randomize = TRUE, multithread=FALSE, nbases = 1e8)
errF <- learnErrors(filtFs, randomize = TRUE, multithread=FALSE, nbases = 1e8)
plotErrors(errF, nominalQ=TRUE)
## dereplication
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
## Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names
## run the core sample inference program
dadaFs <- dada(derepFs, err=errF, multithread=FALSE)
dadaRs <- dada(derepRs, err=errR, multithread=FALSE)
## insepct one
dadaFs[[1]]
## merge paired reads
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])
## construct amplicon sequence variant table (ASV), an OTU table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
## inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))
## remove non-target lengths 
#seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% seq(250,254)]
## table(nchar(getSequences(seqtab2)))
## identify chimers
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=FALSE, verbose=TRUE)
dim(seqtab.nochim)
## define function getN
getN <- function(x) sum(getUniques(x))
## track reads through pipeline
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
## see how much data you have left
head(track)
trackRun <- as.matrix(colSums(track))
## calculate Jaccard coefficients
Rich <- ncol(seqtab.nochim)
N <- seqtab.nochim
N[N>0] <- 1
N[is.na(N)] <- 0
k <- philentropy::distance(N, method = "jaccard")
# write.csv(k, "k.csv", row.names = rownames(N))
G12 <- 1-k[1,4]
H12 <- 1-k[2,3]
GH <- mean(G12, H12)
NoChim <- sum(seqtab.nochim)/sum(seqtab)
## download taxonmy
dna <- DNAStringSet(getSequences(seqtab.nochim)) # Create a DNAStringSet from the ASVs
load("C:/Users/mglam/Dropbox/BIOINFO/Rscripts/SILVA_SSU_r138_2019.RData") # CHANGE TO THE PATH OF YOUR TRAINING SET
# Match reads against reference sequencies
ids <- IdTaxa(dna, trainingSet, strand="top", processors=NULL, verbose=TRUE, threshold = 0) # use all processors
ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species") # ranks of interest
# Convert the output object of class "Taxa" to a matrix analogous to the output from assignTaxonomy
taxid <- t(sapply(ids, function(x) {
  m <- match(ranks, x$rank)
  taxa <- x$taxon[m]
  taxa[startsWith(taxa, "unclassified_")] <- NA
  taxa
}))
colnames(taxid) <- ranks; rownames(taxid) <- getSequences(seqtab.nochim)
# Extract matrix of confidence intervals 
output <- t(sapply(ids, function(x) {
  m <- match(ranks, x$rank)
  confd <- x$confidence[m]
  confd
}))
# Calculate average bootstrap value for domain - class
BootClass <- sum(output[,1:3], na.rm=TRUE)/(3*Rich) 
BootGenus <- sum(output[,1:6], na.rm=TRUE)/(6*Rich) 
## name each data value
SNRalgn <- c(try, t1, t2, t3, Rich, G12, H12, GH, NoChim, M12, BootClass, BootGenus)
names(SNRalgn) <- c("try", "t1", "t2", "t3", "Rich", "G12", "H12", "Rep", "NoChim", "M12", "BootClass", "BootGenus")
B = matrix(SNRalgn, nrow=12, ncol=1) 
# add number of reads at processing
E <- rbind(B,trackRun)
D <- cbind(E, D)
rownames(D) <- c("try", "t1", "t2", "t3", "Rich", "G12", "H12", "GH", "NoChim", "M12", "BootClass", "BootGenus", "input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim" )
M12 <- max(D[8,])
#saveRDS(D, "dadaD2optmz.rds")
plot.default(D[1,], D[10,], xlab = "Try", ylab = "Duplicate")
if (length(D) > 200) {break}
}
```

    ## 495690 total bases in 3813 reads from 4 samples will be used for learning the error rates.
    ## 713031 total bases in 3813 reads from 4 samples will be used for learning the error rates.
    ## Sample 1 - 951 reads in 550 unique sequences.
    ## Sample 2 - 962 reads in 558 unique sequences.
    ## Sample 3 - 956 reads in 517 unique sequences.
    ## Sample 4 - 944 reads in 578 unique sequences.
    ## Sample 1 - 951 reads in 376 unique sequences.
    ## Sample 2 - 962 reads in 451 unique sequences.
    ## Sample 3 - 956 reads in 446 unique sequences.
    ## Sample 4 - 944 reads in 479 unique sequences.
    ## ================================================================================
    ## 
    ## Time difference of 12.39 secs

![](dada2optimize_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

    ## 687050 total bases in 3775 reads from 4 samples will be used for learning the error rates.
    ## 702150 total bases in 3775 reads from 4 samples will be used for learning the error rates.
    ## Sample 1 - 947 reads in 545 unique sequences.
    ## Sample 2 - 950 reads in 542 unique sequences.
    ## Sample 3 - 949 reads in 510 unique sequences.
    ## Sample 4 - 929 reads in 563 unique sequences.
    ## Sample 1 - 947 reads in 478 unique sequences.
    ## Sample 2 - 950 reads in 537 unique sequences.
    ## Sample 3 - 949 reads in 522 unique sequences.
    ## Sample 4 - 929 reads in 595 unique sequences.
    ## ================================================================================
    ## 
    ## Time difference of 15.65 secs

![](dada2optimize_files/figure-gfm/unnamed-chunk-8-2.png)<!-- -->

    ## 465309 total bases in 3783 reads from 4 samples will be used for learning the error rates.
    ## 722553 total bases in 3783 reads from 4 samples will be used for learning the error rates.
    ## Sample 1 - 942 reads in 550 unique sequences.
    ## Sample 2 - 954 reads in 551 unique sequences.
    ## Sample 3 - 953 reads in 518 unique sequences.
    ## Sample 4 - 934 reads in 577 unique sequences.
    ## Sample 1 - 942 reads in 357 unique sequences.
    ## Sample 2 - 954 reads in 433 unique sequences.
    ## Sample 3 - 953 reads in 437 unique sequences.
    ## Sample 4 - 934 reads in 448 unique sequences.
    ## ================================================================================
    ## 
    ## Time difference of 14.2 secs

![](dada2optimize_files/figure-gfm/unnamed-chunk-8-3.png)<!-- -->

    ## 611200 total bases in 3820 reads from 4 samples will be used for learning the error rates.
    ## 676140 total bases in 3820 reads from 4 samples will be used for learning the error rates.
    ## Sample 1 - 955 reads in 539 unique sequences.
    ## Sample 2 - 959 reads in 539 unique sequences.
    ## Sample 3 - 958 reads in 504 unique sequences.
    ## Sample 4 - 948 reads in 561 unique sequences.
    ## Sample 1 - 955 reads in 463 unique sequences.
    ## Sample 2 - 959 reads in 527 unique sequences.
    ## Sample 3 - 958 reads in 508 unique sequences.
    ## Sample 4 - 948 reads in 579 unique sequences.
    ## ================================================================================
    ## 
    ## Time difference of 15.04 secs

![](dada2optimize_files/figure-gfm/unnamed-chunk-8-4.png)<!-- -->

    ## 563270 total bases in 3634 reads from 4 samples will be used for learning the error rates.
    ## 770408 total bases in 3634 reads from 4 samples will be used for learning the error rates.
    ## Sample 1 - 914 reads in 540 unique sequences.
    ## Sample 2 - 923 reads in 548 unique sequences.
    ## Sample 3 - 918 reads in 510 unique sequences.
    ## Sample 4 - 879 reads in 551 unique sequences.
    ## Sample 1 - 914 reads in 427 unique sequences.
    ## Sample 2 - 923 reads in 474 unique sequences.
    ## Sample 3 - 918 reads in 461 unique sequences.
    ## Sample 4 - 879 reads in 516 unique sequences.
    ## ================================================================================
    ## 
    ## Time difference of 15.37 secs

![](dada2optimize_files/figure-gfm/unnamed-chunk-8-5.png)<!-- -->

    ## 589310 total bases in 3802 reads from 4 samples will be used for learning the error rates.
    ## 707172 total bases in 3802 reads from 4 samples will be used for learning the error rates.
    ## Sample 1 - 951 reads in 549 unique sequences.
    ## Sample 2 - 958 reads in 554 unique sequences.
    ## Sample 3 - 953 reads in 512 unique sequences.
    ## Sample 4 - 940 reads in 575 unique sequences.
    ## Sample 1 - 951 reads in 452 unique sequences.
    ## Sample 2 - 958 reads in 504 unique sequences.
    ## Sample 3 - 953 reads in 491 unique sequences.
    ## Sample 4 - 940 reads in 570 unique sequences.
    ## ================================================================================
    ## 
    ## Time difference of 16.34 secs

![](dada2optimize_files/figure-gfm/unnamed-chunk-8-6.png)<!-- -->

    ## 673839 total bases in 3807 reads from 4 samples will be used for learning the error rates.
    ## 658611 total bases in 3807 reads from 4 samples will be used for learning the error rates.
    ## Sample 1 - 953 reads in 530 unique sequences.
    ## Sample 2 - 957 reads in 522 unique sequences.
    ## Sample 3 - 956 reads in 490 unique sequences.
    ## Sample 4 - 941 reads in 546 unique sequences.
    ## Sample 1 - 953 reads in 477 unique sequences.
    ## Sample 2 - 957 reads in 540 unique sequences.
    ## Sample 3 - 956 reads in 525 unique sequences.
    ## Sample 4 - 941 reads in 600 unique sequences.
    ## ================================================================================
    ## 
    ## Time difference of 23.77 secs

![](dada2optimize_files/figure-gfm/unnamed-chunk-8-7.png)<!-- -->

    ## 588896 total bases in 3824 reads from 4 samples will be used for learning the error rates.
    ## 676848 total bases in 3824 reads from 4 samples will be used for learning the error rates.
    ## Sample 1 - 955 reads in 539 unique sequences.
    ## Sample 2 - 962 reads in 542 unique sequences.
    ## Sample 3 - 958 reads in 505 unique sequences.
    ## Sample 4 - 949 reads in 563 unique sequences.
    ## Sample 1 - 955 reads in 455 unique sequences.
    ## Sample 2 - 962 reads in 506 unique sequences.
    ## Sample 3 - 958 reads in 493 unique sequences.
    ## Sample 4 - 949 reads in 574 unique sequences.
    ## ================================================================================
    ## 
    ## Time difference of 15.42 secs

![](dada2optimize_files/figure-gfm/unnamed-chunk-8-8.png)<!-- -->

    ## 460266 total bases in 3742 reads from 4 samples will be used for learning the error rates.
    ## 748400 total bases in 3742 reads from 4 samples will be used for learning the error rates.
    ## Sample 1 - 932 reads in 543 unique sequences.
    ## Sample 2 - 946 reads in 553 unique sequences.
    ## Sample 3 - 947 reads in 518 unique sequences.
    ## Sample 4 - 917 reads in 570 unique sequences.
    ## Sample 1 - 932 reads in 353 unique sequences.
    ## Sample 2 - 946 reads in 426 unique sequences.
    ## Sample 3 - 947 reads in 430 unique sequences.
    ## Sample 4 - 917 reads in 440 unique sequences.
    ## ================================================================================
    ## 
    ## Time difference of 13.27 secs

![](dada2optimize_files/figure-gfm/unnamed-chunk-8-9.png)<!-- -->

    ## 571840 total bases in 3574 reads from 4 samples will be used for learning the error rates.
    ## 793428 total bases in 3574 reads from 4 samples will be used for learning the error rates.
    ## Sample 1 - 900 reads in 539 unique sequences.
    ## Sample 2 - 910 reads in 551 unique sequences.
    ## Sample 3 - 906 reads in 506 unique sequences.
    ## Sample 4 - 858 reads in 539 unique sequences.
    ## Sample 1 - 900 reads in 422 unique sequences.
    ## Sample 2 - 910 reads in 473 unique sequences.
    ## Sample 3 - 906 reads in 457 unique sequences.
    ## Sample 4 - 858 reads in 509 unique sequences.
    ## ================================================================================
    ## 
    ## Time difference of 15.72 secs

![](dada2optimize_files/figure-gfm/unnamed-chunk-8-10.png)<!-- -->

    ## 493739 total bases in 3769 reads from 4 samples will be used for learning the error rates.
    ## 738724 total bases in 3769 reads from 4 samples will be used for learning the error rates.
    ## Sample 1 - 936 reads in 545 unique sequences.
    ## Sample 2 - 952 reads in 551 unique sequences.
    ## Sample 3 - 950 reads in 518 unique sequences.
    ## Sample 4 - 931 reads in 576 unique sequences.
    ## Sample 1 - 936 reads in 368 unique sequences.
    ## Sample 2 - 952 reads in 442 unique sequences.
    ## Sample 3 - 950 reads in 442 unique sequences.
    ## Sample 4 - 931 reads in 471 unique sequences.
    ## ================================================================================
    ## 
    ## Time difference of 12.57 secs

![](dada2optimize_files/figure-gfm/unnamed-chunk-8-11.png)<!-- --> \#\#
Export data to data frame and plot results

``` r
Dd <- as.data.frame(t(D))
plot.default(Dd$try, Dd$M12, xlab = "try", ylab = "Max Duplicate")
```

![](dada2optimize_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

## Export data to data frame and plot results

``` r
Dc <- readRDS("De.rds")
De <- rbind(Dd, Dc)
# saveRDS(De, "De.rds")
```

## Plot pairs

``` r
n <- c(5,8,12,18)
pairs(De[,n], pch = 19)
```

![](dada2optimize_files/figure-gfm/unnamed-chunk-11-1.png)<!-- --> \#\#
Plot number of ASVs versus duplication rate

``` r
plot.default(De$Rich, De$GH, xlab = "ASV", ylab = "Duplicate")
```

![](dada2optimize_files/figure-gfm/unnamed-chunk-12-1.png)<!-- --> \#\#
Plot number of ASVs versus duplication rate

``` r
plot.default(De$Rich, De$BootGenus, xlab = "ASV", ylab = "BootGenus")
```

![](dada2optimize_files/figure-gfm/unnamed-chunk-13-1.png)<!-- --> \#\#
Plot duplication versus match to confidence of assignment to genus

``` r
plot.default(De$BootGenus, De$GH, xlab = "Bootstrap values", ylab = "Duplication")
```

![](dada2optimize_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

## Export data to data frame and plot results

``` r
write.csv(as.matrix(t(De)), "De.csv")
plot.default(De$nonchim, De$BootGenus, xlab = "reads", ylab = "Confidence")
```

![](dada2optimize_files/figure-gfm/unnamed-chunk-15-1.png)<!-- --> \#\#
Select parameters

``` r
DeFit <- De[ which(De$Rich > 100),] # set ASV mininum
DeFit <- DeFit[order(-DeFit$GH, -DeFit$BootGenus),] # sort by fit
t1 <- DeFit$t1[1]
t2 <- DeFit$t2[1]
t3 <- DeFit$t3[1]
write.csv(as.matrix(t(DeFit)), "dadaFit.csv")
tDF <- c(DeFit$t3[1], DeFit$t1[1], DeFit$t2[1])
names(tDF) <- c("trimLeft", "truncateLeft", "truncateRight")
tDF = as.matrix(tDF, nrow=3, ncol=1) 
colnames(tDF) <- c("parameter")
tDF
```

    ##               parameter
    ## trimLeft             18
    ## truncateLeft        220
    ## truncateRight       162

End
