# Day 4 - Data Analysis, Part 1: The basics of Differential Expression analysis with DESeq2
## Friday, November 6th, 2020


<br>

### Join the Slack channels! 
#### https://join.slack.com/t/cvrcbioinform-zmq7258/shared_invite/zt-hja1edwl-RPS1PxgXixBjfUStcW_JTg

<br>  

##

This session will focus on basic RNAseq analysis using the DESeq2 R package, using the same sample RNAseq dataset we used in the past two sessions (Data Exploration [Part 1](https://github.com/florschlamp/CVRC_NYU_Langone/tree/master/Bioinformatics_Workshop_Series_2020/Day_2) and [Part 2](https://github.com/florschlamp/CVRC_NYU_Langone/tree/master/Bioinformatics_Workshop_Series_2020/Day_3)). First we will show you how to do a quick data exploration, including plotting PCAs, with built in DESeq2 functions. Next, we will walk you through performing Differential Expression analysis, how to get a results table, and explore it using p-value and log Fold Change cutoffs. You will learn to extract the most significantly expressed genes and plot them based on their log Fold Change. You will also learn to extract the normalized counts from the DESeq2 package, and export them along with the DE results as files you can later access with Excel.

<br>

### REQUIREMENTS TO PARTICIPATE IN THIS SESSION:
 
Make sure you have the following packages installed in Rstudio: DESeq2, ggplot2.  
Example of how to install R packages:
```r
### Installing DESeq2 from Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")

# if the installation worked, you should be able to run the following line without error messages:
library(DESeq2)


### Installing ggplot2: 
install.packages("ggplot2")

# if the installation worked, you should be able to run the following line without error messages:
library(ggplot2)
```

### OFFICE HOURS:
 
Friday Nov. 6th
* 11am to noon  
* 1:30pm until the start of the session at 2:30pm  
 
Zoom link:
https://nyulangone.zoom.us/j/98456141720
