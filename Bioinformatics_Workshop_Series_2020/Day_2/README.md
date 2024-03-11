# Day 2 - Data Exploration, Part 1: Basic RNA-seq Data Manipulation
## Friday, October 9th, 2020


<br>

### Join the Slack channels! 
#### https://join.slack.com/t/cvrcbioinform-zmq7258/shared_invite/zt-hja1edwl-RPS1PxgXixBjfUStcW_JTg

<br>  

##

In this session we will use a sample RNAseq dataset (raw counts and metadata) to go over the basic steps of RNA-seq data exploration. This includes reading in data from a file and performing initial exploratory analysis such as calculating library size for each sample and plotting them as barplots. You will learn to do basic data processing of the raw counts, such as removing genes that have zero counts across all samples, performing basic normalization and log transformation, and filter genes with low read counts. By the end of this session you will also be able to plot PCAs of the data, and use this information to remove outliers.


<br>

### REQUIREMENTS TO PARTICIPATE IN THIS SESSION:

Make sure you have the following packages installed in Rstudio: ggplot2 and edgeR (note that edgeR needs to be downloaded through Bioconductor). 

```r
# installing ggplot2: 
install.packages("ggplot2")


# installing edgeR:
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("edgeR")


# if the installations worked, you should be able to run the following lines without error messages:
library(ggplot2)
library(edgeR)
```

### OFFICE HOURS:
 
Friday Oct. 9th
* 11am to noon  
* 1:30pm until the start of the session at 2:30pm  


