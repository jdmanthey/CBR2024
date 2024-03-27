## Final Project Template (may be some small errors still)

### 1. Get Ready

If you are using the same computer as always, you should be good to go. If not, you may need to install any 
necessary packages as shown in the previous tutorials.

### 2. Download data

You will want to create a working directory (on your desktop, on a flash drive, all that matters is that you know what it is
called and where it is). In this working directory, create a new directory called: "raw_data" and add the gg_13_8_train_set_97.fa.gz file.
When you have your samples downloaded (just download your subset), you will place all these fastq.gz files in the "raw_data" folder. 

### 3. Load packages one at a time to make sure they work.

Make sure after each "library" command that you do not get any errors otherwise you will run into problems later.

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
    
    library(vegan)
    
    library(dplyr)
    
    library(microbiome)


### 4. Set up analysis

Set your working directory to that which you created. 

Use list.files() to check that you have the directory raw_data in your current path as well as the 
gg_13_8_train_set_97.fa.gz file. If not, your path is incorrect and you need to reset the working directory.

    # set the path to the location of the sequencing reads
    path <- "raw_data"
    
    list.files(path)
    
The files should have been listed with the list.files(path) command. If not, the directories are incorrect.

    # sort the order of the forward and reverse reads
    fnFs <- sort(list.files(path, pattern="_R1.fastq.gz"))

    fnRs <- sort(list.files(path, pattern="_R2.fastq.gz"))

    # extract sample names from files
    sample.names <- sapply(strsplit(fnFs, "_"), `[`, 1)
    
    sample.names

You should be able to see the sample names from the second command above.

    # Specify the full path to the fnFs and fnRs
    fnFs <- file.path(path, fnFs)
    
    fnRs <- file.path(path, fnRs)
    
### 5. Error profiling

    # Look at a few sequences and check out their error profiles
    plotQualityProfile(fnFs[1:3])
    
    plotQualityProfile(fnRs[1:3])
    
    # file paths for putting the filtered reads
    filt_path <- file.path(path, "filtered")
    
    filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
    
    filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))

    # filter and trim the samples
    out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(140,140),
              maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
              compress=TRUE, multithread=FALSE)
              
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
    
    unname(head(taxa))

Save progress:

    save.image(file="microbe_workflow2.Rdata")

### 10. Create phylogenetic tree of microbial data

    seqs <- getSequences(seqtab.nochim)
    
    names(seqs) <- seqs
    
    alignment <- AlignSeqs(DNAStringSet(seqs))
    
    phang.align <- phyDat(as(alignment, "matrix"), type="DNA")
    
    dm <- dist.ml(phang.align)
    
    treeNJ <- NJ(dm)
    
    plot(treeNJ, show.tip.label=F)

Save progress:

    save.image(file="microbe_workflow3.Rdata")

### 11. Make a phyloseq object with the sample data, phylogeny, and sequence variant table

    # Make a data.frame holding the sample data
    samples.out <- rownames(seqtab.nochim)
    
You will have to make a custom data frame based on the sequences you have. For each characteristic of your data you are 
interested in, you should make a vector for each item first, then combine them all into a dataframe. An example is shown here 
with 10 samples, the sample IDs, and three additional traits. You should modify these, rename the traits,
add to or reduce vector size lengths as necessary given your sample size and the traits of interest
    
    sample.id <- sample.names
    
    trait1 <- c("a", "b", "c", "d", "e", "f", "g", "h", "i", "j")
    
    trait2 <- c(1, 3, 1, 1, 5, 3, 5, 4, 4, 5)
    
    trait3 <- c("a", "a", "d", "d", "f", "f", "g", "g", "j", "j")
    
    micro.df <- data.frame(Sample.ID=sample.id, trait_name1=trait1, trait_name2=trait2, trait_name3=trait3)

    rownames(micro.df) <- samples.out


    # construct a phyloseq object
    ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=F),
    		sample_data(micro.df),
    		tax_table(taxa),
    		phy_tree(treeNJ))
        
    ps
    
    # remove cyanobacteria
    ps <- subset_taxa(ps, Phylum != "p__Cyanobacteria")

Save progress:

    save.image(file="microbe_workflow4.Rdata")

### 12. Composition

Look at a summary of the numbers of unique sequence variants:
    
    summary(sample_sums(ps))

Find the unique phyla, classes, and families in the dataset:

    sort(get_taxa_unique(ps, "Phylum"))
    
    sort(get_taxa_unique(ps, "Class"))
    
    sort(get_taxa_unique(ps, "Family"))

The following lines of code will summarize the phyla and classes per sample, and summaries across all samples, followed by
writing those results to two tables each (copy and paste the whole thing):

    i_phylum <- c()
    reps <- as.vector(unique(ps@tax_table[ ,colnames(ps@tax_table) == "Phylum"]))
    reps_matrix <- ps@tax_table[,colnames(ps@tax_table) == "Phylum"]
    for(a in 1:length(reps)) {
    	a.rep <- rownames(reps_matrix)[reps_matrix == reps[[a]]]
    	a.rep <- apply(ps@otu_table[ ,match(a.rep, colnames(ps@otu_table))], 1, sum)
    	i_phylum <- cbind(i_phylum, a.rep)
    }
    colnames(i_phylum) <- reps
    inds_names <- rownames(i_phylum)
    i_class <- c()
    reps <- as.vector(unique(ps@tax_table[ ,colnames(ps@tax_table) == "Class"]))
    reps <- na.omit(reps)
    reps_matrix <- ps@tax_table[,colnames(ps@tax_table) == "Class"]
    reps_matrix <- na.omit(reps_matrix)
    for(a in 1:length(reps)) {
    	a.rep <- rownames(reps_matrix)[reps_matrix == reps[[a]]]
    	a.rep <- apply(ps@otu_table[ ,match(a.rep, colnames(ps@otu_table))], 1, sum)
    	i_class <- cbind(i_class, a.rep)
    }
    colnames(i_class) <- reps

    # phylum output					 
    i_tmp <- c()
    for(a in 1:nrow(i_phylum)) {
    	a.rep <- i_phylum[a,]
    	a.rep <- a.rep / sum(a.rep) * 100
    	i_tmp <- rbind(i_tmp, a.rep)
    }
    i_phylum <- i_tmp
    i_phylum_output <- cbind(colnames(i_phylum), apply(i_phylum, 2, mean), apply(i_phylum, 2, min), apply(i_phylum, 2, max))
    i_phylum_output <- data.frame(Phylum=as.vector(sapply(strsplit(i_phylum_output[,1], "__"), "[[", 2)), Mean=(as.numeric(i_phylum_output[,2])),Min=(as.numeric(i_phylum_output[,3])), Max=(as.numeric(i_phylum_output[,4])))

    # class output
    i_tmp <- c()
    for(a in 1:nrow(i_class)) {
    	a.rep <- i_class[a,]
    	a.rep <- a.rep / sum(a.rep) * 100
    	i_tmp <- rbind(i_tmp, a.rep)
    }
    i_class <- i_tmp
    i_class_output <- cbind(colnames(i_class), apply(i_class, 2, mean), apply(i_class, 2, min), apply(i_class, 2, max))
    i_class_output <- data.frame(Class=substr(i_class_output[,1], 4, nchar(i_class_output[,1])), Mean=(as.numeric(i_class_output[,2])),Min=(as.numeric(i_class_output[,3])), Max=(as.numeric(i_class_output[,4])))

    write.table(i_phylum_output, file="i_phylum_output.txt", sep="\t", quote=F, row.names=F)
    write.table(i_class_output, file="i_class_output.txt", sep="\t", quote=F, row.names=F)
    write.table(cbind(inds_names, i_phylum), file="i_phylum_output2.txt", sep="\t", quote=F, row.names=F)
    write.table(cbind(inds_names, i_class), file="i_class_output2.txt", sep="\t", quote=F, row.names=F)

You can look at what the summary tables look like here:
  
    i_phylum_output
    i_class_output

Now we can plot the phyla that are found in the samples:

    sample_number <- as.data.frame(sample_data(ps)$Sample.ID)
    colnames(sample_number) <- c("sample_number")
    phyla <- melt(cbind(sample_number, as.data.frame(i_phylum)), id=c("sample_number"), value.name="phyla")
    ggplot(data=phyla, aes(x=sample_number, y=value, fill=variable)) + geom_bar(stat="identity") + scale_fill_manual(values=(turbo(length(unique(phyla$variable)))))

Do the same plot, but for classes:

    class <- melt(cbind(sample_number, as.data.frame(i_class)), id=c("sample_number"), value.name="phyla")
    ggplot(data=class, aes(x=sample_number, y=value, fill=variable)) + geom_bar(stat="identity") + scale_fill_manual(values=(turbo(length(unique(class$variable)))))

### 13. Alpha diversity

Make a table of many types of alpha diversity:

    alpha_table <- estimate_richness(ps)
    
    alpha_table

Let's make a simple table that includes 4 measures of alpha diversity. We'll save them to the object named 'alpha.'
	
    alpha <- cbind(estimate_richness(ps)[,c(1,6,7)], pd((otu_table(ps)@.Data), ps@phy_tree, include.root=F)$PD)
  
    colnames(alpha) <- c("Observed_SVs", "Shannon", "Simpson", "PD")
  
    alpha

Remember that sampling depth can influence estimates of diversity, especially if they are incompletely sampled communities.
Let's plot some rarefaction curves to take one look at this factor:

	par(mar=c(4.5,4.5,3,3)) # plotting margins
	rarecurve_table <- otu_table(ps)
  	class(rarecurve_table) <- "matrix"
  	rarecurve(rarecurve_table, step=50, cex=0.5, label=F, xlim=c(0,50000))
  	
Let's plot a couple types of alpha diversity by trait_name3 (you can use the name of any of your traits here (if they are 
categorical):

    plot_richness(ps, x="trait_name3", measures=c("Shannon", "Simpson"))
     
We can also plot the relationships among each of the estimates of alpha diversity. Here, Observed_SVs = number of observed
sequence variants and PD = phylogenetic diversity.

    plot(alpha, pch=19, cex=1)

### 14. Beta diversity

Now we'll measure beta diversity among communities in a few different ways, colored by trait3 (again, any of your traits will
work. 

First, we'll make a plot of the Bray-Curtis distance:
	
    ord.nmds.bray <- ordinate(ps, method="NMDS", distance="bray")

    plot_ordination(ps, ord.nmds.bray, color="trait_name3", shape="trait3", title="Bray NMDS") + geom_point(size=3)

And unweighted Unifrac distance:

    ord.unweighted.unifrac <- ordinate(ps, "PCoA", "unifrac", weighted=F)

    plot_ordination(ps, ord.unweighted.unifrac, color="trait_name3", shape="trait3",title="Unweighted Unifrac PCoA") + geom_point(size=3)
	
And weighted Unifrac distance:

    ord.weighted.unifrac <- ordinate(ps, "PCoA", "unifrac", weighted=T)

    plot_ordination(ps, ord.weighted.unifrac, color="trait_name3", shape="trait3",title="Weighted Unifrac PCoA") + geom_point(size=3)

Save progress:

    save.image(file="microbe_workflow5.Rdata")

### 15. Stats

Use a relevant statistical test that was covered in class (or if you want, one we haven't covered) to test your hypotheses.
This will be different for each dataset and question, and will likely involve a little bit of trouble shooting. Feel free
to ask questions or search the internet for solutions. I expect this may take up more time than many of the other parts
of the project.

### 16. Subset a phyloseq object to the core microbiome

May be useful for some people's projects (you choose the prevalence (or different ranges of prevalences) that you want):

    ps_core <- core(ps, detection = 0, prevalence = .3) 
    


    

