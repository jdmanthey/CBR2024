## Some statistics with microbiome data 29 February 2024

Today we will use R to perform some statistical analyses on microbiome datasets. This should give you an idea about some 
of the possibilities for hypothesis testing given the data types available.

I'd suggest taking a look at lecture notes at the same time as performing analyses if you think it will help with
understanding and retention.

### 1. Set up working directory and download data

Go to the same working directory used for the other microbiome analyses that you have done in class so far. You will need to
load your Rhistory if it doesn't automatically populate, or you can load the object called "microbe_workflow_day2.Rdata" 
If you do not have the saved data for whatever reason, I have made a copy that I used myself. Download it 
[here](https://drive.google.com/file/d/1E0_NrX0AsbOpX71Qek9tYMRXy0TJdUCz/view?usp=drive_link) 
and put it in your working directory.

You will also need to download a small table of host genetic diversity for one of the exercises. Download that table 
[here](https://drive.google.com/open?id=1sXynwQm-UAsuM0G6sKuGLUv42763vkCb) and put it in your working directory.

### 2. Load packages one at a time to make sure they work

    library(phyloseq)

    library(ggplot2)

    library(vegan)

### 3. T-test

We'll do a T-test of phylogenetic diversity between the Chiricahua and Huachuca populations of hosts.

First, subset the phylogenetic diversity into two vectors:

    chiricahua_pd <- alpha$PD[sample_data(ps)$Location == "Chiricahua"]

    huachuca_pd <- alpha$PD[sample_data(ps)$Location == "Huachuca"]

Now we can run the t-test:

    t.test(chiricahua_pd, huachuca_pd)

### 4. ANOVA

Now we will run an ANOVA with the phylogenetic diversity for all populations.

First, we will set up our data frame with two columns and take a look at it:

    all_pd <- data.frame(PD=as.numeric(alpha$PD), population=sample_data(ps)$Location)

    all_pd
    
Now we will run the ANOVA, trying to explain the phylogenetic diversity by population information. Note the structure of the 
equation in the aov code (response_variable ~ explanatory_variable).

    pd_aov_out <- aov(PD ~ population, data=all_pd)

And summarize the anova model output.

    summary(pd_aov_out)
    
And finally run a Tukey Test of the ANOVA output:

    TukeyHSD(pd_aov_out)

### 5. Linear Regression

We will use linear regression to look at the relationship between host genetic diversity (measured by observed heterozygosity)
and alpha diversity of their microbial communities.

First, read in the host diversity table (that you downloaded in step #1):

    host_diversity <- read.table("host_diversity.txt", sep="\t", header=T, stringsAsFactors=F)

Plot host diversity against the microbial community phylogenetic diversity. Do you expect a significant relationship between the two variables?

    plot(host_diversity$heterozygosity, all_pd$PD, pch=19, xlab="Host Heterozygosity", ylab="Microbiome Phylogenetic Diversity")

Run a linear model of the two variables:

    lm_out <- lm(all_pd$PD ~ host_diversity$heterozygosity)

Summarize the output:

    summary(lm_out)
    
You can also look at the predicted relationship by adding a line to your plot:

    abline(lm_out)

### 6. ADONIS test

Let's run an ADONIS test to look at the relationship between differences in OTU abundance in the microbial communities as a
function of the location and genetic diversity of the hosts:

    adonis2(otu_table(ps) ~ sample_data(ps)$Location + host_diversity$heterozygosity)

