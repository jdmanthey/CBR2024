## Core Microbiome 19 March 2024

Today we will use R to identify the core microbiome in our test dataset. 

### 1. Set up working directory and download data

Go to the same working directory used for the other microbiome analyses that you have done in class so far. You will need to
load the object called "microbe_workflow_day2.Rdata." If you do not have the saved data for whatever reason, I have made a 
copy that I used myself. Download it [here](https://drive.google.com/file/d/1E0_NrX0AsbOpX71Qek9tYMRXy0TJdUCz/view?usp=drive_link) and put it 
in your working directory.

    load(file="microbe_workflow_day2.Rdata")

### 2. Load packages one at a time to make sure they work

    library(phyloseq)

    library(ggplot2)

    library(vegan)
    
    library(phyloseq)

    library(ggplot2)

    library(RColorBrewer)

    library(vegan)

Install and load the new necessary packages (one line at a time):

    install.packages("dplyr")

    library(dplyr)
    
    install.packages("devtools")

    library("devtools")

    install_github("microbiome/microbiome")

    library(microbiome)


### 3. Identify the core microbiome

Transform the phyloseq counts data to relative abundance:

    ps_rel <- microbiome::transform(ps, "compositional")

Rename the sequence variant names with numbers and taxonomy identifiers:

    temp_taxa_names <- paste(as.vector(tax_table(ps_rel)[match(taxa_names(ps_rel), rownames(tax_table(ps_rel))),5]), ":", as.vector(tax_table(ps_rel)[match(taxa_names(ps_rel), rownames(tax_table(ps_rel))),6]), sep="")
    
    temp_taxa_names <- paste(seq(from=1, to=length(temp_taxa_names)), ":", temp_taxa_names, sep="")
    
    taxa_names(ps_rel) <- temp_taxa_names
    
What are the 6 most prevalent bacteria?

    head(prevalence(ps_rel, detection = 1/100, sort = TRUE))

What are the core members of the microbial communities at 30% prevalence across samples?

    # Command to rerun for question 4 on worksheet
    core_members(ps_rel, detection = 0, prevalence = .3) 

Make a new phyloseq object with only the core microbiota:

    # Command to rerun for question 4 on worksheet
    pseq_core <- core(ps_rel, detection = 0, prevalence = .3) 
    
What is the relative abundance of the core in each of the samples?

    # Command to rerun for question 4 on worksheet
    sample_sums(core(ps_rel, detection = 0, prevalence = .3)) 

### 4. Plot the core microbiome

Define the detection and prevalence limits of the plot:
   
    det <- c(0.5,1, 2, 5, 20)/100

    prevalences <- seq(.05, 1, .05)

And plot:

    # Command to rerun for question 4 on worksheet
    plot_core(pseq_core, prevalences = prevalences, detections = det, plot.type = "heatmap", horizontal=F,colours = rev(brewer.pal(5, "Spectral"))) + xlab("Relative Abundance (%)")


