## Microbiome Community Data Processing - February 2024

For this practical exercise, we will use R with a small practice dataset to learn about filtering and processing
raw microbial community sequencing data. The data are the V4 region of the 16S rRNA gene and come from intestinal microbe
samples of a bird species: _Sitta carolinensis_. The bird gut samples we will use here were sampled from 4 mountain ranges
in southeastern Arizona. If you are interested in looking them up, their names are: (1) Chiricahua, (2) Huachuca, 
(3) Pinaleño, and (4) Santa Catalina Mountain Ranges.

Throughout, you may want to refer to lecture notes about microbiome processing, particularly when looking at some of the 
figures.

### 1. Download data

Download the zip file from this link: [link](https://drive.google.com/file/d/1FXaJmU4JY4f9ibGp4hykT31W7NcN5OS_/view?usp=drive_link), unzip it,
and set that directory (class_subsample, not raw_data) as your working directory in RStudio.

### 2. Install necessary packages

    # We are using the DADA2 package (https://benjjneb.github.io/dada2/)
    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
    BiocManager::install("dada2")

    # We need the DECIPHER package for some functions
    BiocManager::install("DECIPHER")
    
    # And phyloseq for diversity and plotting of the community data
    BiocManager::install("phyloseq")
    
    # We need a couple other packages that are not included in the above installs
    install.packages("reshape")
    install.packages("picante")
    install.packages("viridis")
    install.packages("phangorn")


### 3. Load packages one at a time to make sure they work

    library(dada2)
    library(phangorn)
    library(DECIPHER)
    library(phyloseq)
    library(ggplot2)
    library(reshape)
    library(RColorBrewer)
    library(plyr)
    library(picante)
    library(viridis)

### 4. Set up of analysis

The data we have today is a subset of the full data collected from the birds. We have 5000 paired-end Illumina sequencing
reads for each individual. These reads target one of the variable regions in the 16S rRNA gene in microorganisms. 

    # set the path to the location of the sequencing reads
    path <- "raw_data"
    list.files(path)
    
The files should have been listed with the second command. If not, the directories are incorrect.

    # sort the order of the forward and reverse reads
    fnFs <- sort(list.files(path, pattern="_R1_001.fastq"))
    fnRs <- sort(list.files(path, pattern="_R2_001.fastq"))

    # extract sample names from files
    sample.names <- paste(sapply(strsplit(fnFs, "_"), `[`, 1), sapply(strsplit(fnFs, "_"), `[`, 4), sapply(strsplit(fnFs, "_"), `[`, 5), sep="_")
    sample.names

You can see we have file names for 11 individuals and that these names include the locality information.

    # Specify the full path to the fnFs and fnRs
    fnFs <- file.path(path, fnFs)
    fnRs <- file.path(path, fnRs)
    
### 5. Error profiling

    # Look at a few sequences and check out their error profiles
    plotQualityProfile(fnFs[1:2])
    plotQualityProfile(fnRs[1:2])
    
    # file paths for putting the filtered reads
    filt_path <- file.path(path, "filtered")
    filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
    filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))

    # filter and trim the samples
    out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(140,140),
              maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
              compress=TRUE, multithread=FALSE)
              
    # look at the top few lines of the output of the filtering and trimming
    head(out)
    
    # learn the error rates for the sequencing
    errF <- learnErrors(filtFs, multithread=FALSE)
    errR <- learnErrors(filtRs, multithread=FALSE)

    # look at plots of errors
    plotErrors(errF, nominalQ=TRUE)
    plotErrors(errR, nominalQ=TRUE)

### 6. Dereplicate and call sequence variants

    # dereplicate all reads
    derepFs <- derepFastq(filtFs, verbose=TRUE)
    derepRs <- derepFastq(filtRs, verbose=TRUE)
    
    # rename the dereplicated reads files
    names(derepFs) <- sample.names
    names(derepRs) <- sample.names

    # run the main file to call all of the unique sequences
    dadaFs <- dada(derepFs, err=errF, pool=T,multithread=TRUE)
    dadaRs <- dada(derepRs, err=errR, pool=T,multithread=TRUE)

Save progress:

    save.image(file="microbe_workflow1.Rdata")

### 7. Merge forward/reverse reads, and remove chimeras

    # merge the forward and reverse sequences
    mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
    
    # make a table of all sequences
    seqtab <- makeSequenceTable(mergers)
    
    # Inspect distribution of sequence lengths
    table(nchar(getSequences(seqtab)))
    
    # keep all mergers with length near the mode (253)
    seqtab <- seqtab[,nchar(colnames(seqtab)) %in% seq(252,254)]
    
    # remove chimeras
    seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
    
    # proportion of sequences not chimeric
    sum(seqtab.nochim)/sum(seqtab)

### 8. Summarize filtering

    # summarize the filtering
    getN <- function(x) sum(getUniques(x))

    # make a new matrix with all the counts at different stages
    track <- cbind(out, sapply(dadaFs, getN), sapply(mergers, getN), rowSums(seqtab), rowSums(seqtab.nochim))

    # add column and row names to the matrix
    colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled", "nonchim")
    rownames(track) <- paste0(sapply(strsplit(sample.names, "_"), "[[", 1), "_", sapply(strsplit(sample.names, "_"), "[[", 3))

    # show the matrix
    track

### 9. Assign taxonomy

Here, we will use the GreenGenes database formatted for DADA2 to classify the 16S sequences we have here. Note that the
classification is only accurate to the family level, and any inferences about genus or species may or may not be accurate.

    taxa <- assignTaxonomy(seqtab.nochim, "gg_13_8_train_set_97.fa.gz", multithread=TRUE)
    
    # what does the output look like? Check out the top few lines:
    unname(head(taxa))

Save progress:

    save.image(file="microbe_workflow1.Rdata")
