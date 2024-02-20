## Microbiome Community Abundance and Diversity 20 February 2024

We will continue from the preliminary filtering of microbiome data we previously did in class. Here, you will
be plotting and interpreting information about the communities' composition and alpha and beta diversity measures.

You should be able to load your Rhistory file. If that doesn't work, in RStudio you can move to your project
directory and enter the command: load("microbe_workflow1.Rdata")

### 1. Load packages one at a time to make sure they work

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

### 2. Check that many of the previously used objects are in your environment

    # check that the items loaded
    ls()

### 3. Create phylogenetic tree of microbial data

    seqs <- getSequences(seqtab.nochim)
    names(seqs) <- seqs
    alignment <- AlignSeqs(DNAStringSet(seqs))
    phang.align <- phyDat(as(alignment, "matrix"), type="DNA")
    dm <- dist.ml(phang.align)
    treeNJ <- NJ(dm)
    plot(treeNJ, show.tip.label=F)

Save progress:

    save.image(file="microbe_workflow1.Rdata")

### 4. Make a phyloseq object with the sample data, phylogeny, and sequence variant table

    # Make a data.frame holding the sample data
    samples.out <- rownames(seqtab.nochim)
    samp.number <- sapply(strsplit(samples.out, "_"), `[`, 1)
    species.id <- sapply(strsplit(samples.out, "_"), `[`, 2)
    location <- sapply(strsplit(samples.out, "_"), `[`, 3)
    species.location <- paste(species.id, location, sep=".")
    micro.df <- data.frame(Sample.number=samp.number, Species.ID=species.id, Location=location, Species.Location=species.location)
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

    save.image(file="microbe_workflow1.Rdata")

### 5. Composition

Look at a summary of the numbers of unique sequence variants:
    
    summary(sample_sums(ps))

Find the unique phyla, classes, and families in the dataset:

    sort(get_taxa_unique(ps, "Phylum"))
    sort(get_taxa_unique(ps, "Class"))
    sort(get_taxa_unique(ps, "Family"))

Sample a random family:

	sample(sort(get_taxa_unique(ps, "Family")), 1)
	
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

    sample_number <- as.numeric(sample_data(ps)$Sample.number)
    phyla <- melt(as.data.frame(cbind(sample_number, i_phylum)), id=c("sample_number"), value.name="phyla")
    phyla$sample_number <- as.character(phyla$sample_number)
    ggplot(data=phyla, aes(x=sample_number, y=value, fill=variable)) + geom_bar(stat="identity") + scale_fill_manual(values=(turbo(length(unique(phyla$variable)))))

Do the same plot, but for classes:

    class <- melt(as.data.frame(cbind(sample_number, i_class)), id=c("sample_number"), value.name="class")
    class$sample_number <- as.character(class$sample_number)
    ggplot(data=class, aes(x=sample_number, y=value, fill=variable)) + geom_bar(stat="identity") + scale_fill_manual(values=(turbo(length(unique(class$variable)))))

### 6. Alpha diversity

Here, we'll take a look at a couple ways you can investigate alpha diversity. There is a base function in the phyloseq 
package to look at many types of alpha diversity:

	estimate_richness(ps)

However, we'll just be looking at a few to keep things simple. We'll save them to the object named 'alpha.'
	
	alpha <- cbind(estimate_richness(ps)[,c(1,6,7)], pd((otu_table(ps)@.Data), ps@phy_tree, include.root=F)$PD)
	colnames(alpha) <- c("Observed_SVs", "Shannon", "Simpson", "PD")
	alpha

Remember that sampling depth can influence estimates of diversity, especially if they are incompletely sampled communities.
Let's plot some rarefactions curves to take one look at this factor:

	par(mar=c(4.5,4.5,3,3)) # plotting margins
	rarecurve_table <- otu_table(ps)
  	class(rarecurve_table) <- "matrix"
  	rarecurve(rarecurve_table, step=50, cex=0.5, label=F, xlim=c(0,10000))


If we had the full datasets, the picture may look more completely sampled. Anyways, there are ways of using the full dataset
and also ways of rarefying prior to estimating alpha diversity and beta diversity. We will not get into those methods in these
activities.

Let's plot a couple types of alpha diversity:

	plot_richness(ps, x="Location", measures=c("Shannon", "Simpson"))

We can also plot the relationships among each of the estimates of alpha diversity. Here, Observed_SVs = number of observed
sequence variants and PD = phylogenetic diversity.

	plot(alpha, pch=19, cex=1)

### 7. Beta diversity

Now we'll measure beta diversity among communities in a few different ways. 

First, we'll make a plot of the Bray-Curtis distance:
	
	ord.nmds.bray <- ordinate(ps, method="NMDS", distance="bray")
	wu.df <- data.frame(
		x = as.vector(ord.nmds.bray$points[,1]),
		y = as.vector(ord.nmds.bray$points[,2]),
		Species.ID = species.id,
		Species.Location = species.location)
	plot_ordination(ps, ord.nmds.bray, color="Species.Location", shape="Species.Location", title="Bray NMDS") + geom_point(size=3)

And unweighted Unifrac distance:

	ord.unweighted.unifrac <- ordinate(ps, "PCoA", "unifrac", weighted=F)
	wu.df <- data.frame(
		x = ord.unweighted.unifrac$vectors[,1],
		y = ord.unweighted.unifrac$vectors[,2],
		Species.ID = species.id,
		Species.Location = species.location)
	plot_ordination(ps, ord.unweighted.unifrac, color="Species.Location", shape="Species.Location",title="Unweighted Unifrac PCoA") + geom_point(size=3)
	
And weighted Unifrac distance:

	ord.weighted.unifrac <- ordinate(ps, "PCoA", "unifrac", weighted=T)
	wu.df <- data.frame(
		x = ord.weighted.unifrac$vectors[,1],
		y = ord.weighted.unifrac$vectors[,2],
		Species.ID = species.id,
		Species.Location = species.location)
	plot_ordination(ps, ord.weighted.unifrac, color="Species.Location", shape="Species.Location",title="Weighted Unifrac PCoA") + geom_point(size=3)


Save progress:

    save.image(file="microbe_workflow_day2.Rdata")
