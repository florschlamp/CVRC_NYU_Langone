# CVRC Bioinformatics Workshop Series 2020

## Basic Programming and Bioinformatics Analyses using R in Rstudio

### Schedule


#### Coming up:
* Friday September 25th  
Day 1: “Introduction: First Steps on using R in Rstudio”

* Friday October 9th  
Day 2: “Data Exploration, Part 1: Basic RNA-seq Data Manipulation”

* Friday October 23rd  
Day 3: “Data Exploration, Part 2: Basic RNA-seq Plots”

* Friday November 6th  
Day 4: “Data Analysis, Part 1: The basics of Differential Expression analysis with DESeq2”

* Friday November 20th  
Day 5: “Data Analysis, Part 2: Publication-Ready RNAseq Plots”

* Friday December 4th  
Day 6: “Data Analysis, Part 3: Downstream Analyses”

* Friday December 18th  
Day 7: “Troubleshooting and Miscellaneous”


##  
<br> 

### Syllabus
Day 1) “Introduction: First Steps on using R in Rstudio” (Friday September 25th)

We assume zero programming knowledge and zero previous use of R or Rstudio. We will guide you to install R and Rstudio on your personal computers. You will learn to use R libraries, how to find them, install them, and call them for future use. You will learn how to navigate Rstudio, create variables, distinguish between different data types, and perform basic functions on these data types. You will also learn to read files with external data into R and how to perform basic functions on this data.
This session is required for all beginner users to be able to follow along future sessions. Intermediate users are encouraged to use this as a general review and refresher, especially those who haven’t used R in a while or are new to Rstudio.
  
<br>
    
Day 2) “Data Exploration, Part 1: Basic RNA-seq Data Manipulation” (Friday October 9th)

We will use a sample RNAseq dataset (raw counts and metadata) to go over the basic steps of RNA-seq data exploration. This includes reading in data from a file and performing initial exploratory analysis such as calculating library size for each sample and plotting them as barplots. You will learn to do basic data processing of the raw counts, such as removing genes that have zero counts across all samples, performing basic normalization and log transformation, and filter genes with low read counts. By the end of this session you will also be able to plot PCAs of the data, and use this information to remove outliers.

<br>

Day 3) “Data Exploration, Part 2: Basic RNA-seq Plots” (Friday October 23rd)

In this session we will focus on plotting the data that we cleaned up and normalized in the previous session (Part 1, Data Manipulation). You will learn to use R to do pairwise scatterplots to compare replicates and samples. You will learn the basics of the pheatmap R package to make heatmaps of the similarity between samples, and we will introduce the basics of ggplot R package to plot individual genes of interest. These include dot plots to represent gene expression per sample and per group, barplots, and boxplots. You will also learn how to plot multiple genes at the same time or genes in paired data and time series.

<br>

Day 4) “Data Analysis, Part 1: The basics of Differential Expression analysis with DESeq2” (Friday November 6th)

This session will focus on basic RNAseq analysis using the DESeq2 R package, using the same sample RNAseq dataset we used in the past two sessions (Data Exploration Parts 1 and 2). First we will show you how to do a quick data exploration, including plotting PCAs, with built in DESeq2 functions. Next, we will walk you through performing Differential Expression analysis, how to get a results table, and explore it using p-value and log Fold Change cutoffs. You will learn to extract the most significantly expressed genes and plot them based on their log Fold Change. You will also learn to extract the normalized counts from the DESeq2 package, and export them along with the DE results as files you can later access with Excel. 

<br>

Day 5) “Data Analysis, Part 2: Publication-Ready RNAseq Plots” (Friday November 20th)

Now that you have normalized counts, log Fold Changes, and p-values for all genes, you will learn to plot your data using multiple types of plots, such as volcano plots, MA plots, and heatmaps. More specifically, in this session you will learn how to customize these plots in ggplot and pheatmap R packages based on your needs, and we will also introduce a few other helpful R packages for more advanced customization. Examples of plot modifications that you will be able to do by the end of this session include: how to color code significantly up- and down-regulated genes, how to label the name of certain genes of interest, how to make an interactive volcano plot, how to add or remove legends, how to subset data for heatmaps, how to label heatmap columns based on specific groups (such as treatment vs control, or cell types), and how to make scaled heatmaps.

<br>

Day 6) “Data Analysis, Part 3: Downstream Analyses” (Friday December 4th)

In this session we will focus on a few downstream analyses that will allow you to gain more insight from your differential expression analysis. First, you will learn to use the biomaRt R package to extract various types of data on your genes of interest. These include translating gene IDs into gene names and vice versa, extracting specific information on genes such as chromosome location, and mapping genes to their homologues in a different organism (for example, mouse to human). You will also learn to subset and extract genes for functional enrichment analysis, both using the R package clusterProfiler as well as other online tools such as Panther, GOrilla.

<br>

Day 7) “Troubleshooting and Miscellaneous” (Friday December 18th)

In this final session for Fall 2020 we will go over some general troubleshooting for some common problems you might face when programming in R. We will also use this session to revisit any particular issues that people might have had in the past sessions.
