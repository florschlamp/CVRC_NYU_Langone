# Day 6 - Data Analysis, Part 3: Downstream Analyses
## Friday, December 4th, 2020


<br>

### Join the Slack channels! 
#### https://join.slack.com/t/cvrcbioinform-zmq7258/shared_invite/zt-hja1edwl-RPS1PxgXixBjfUStcW_JTg

<br>  

##

In this session we will focus on a few downstream analyses that will allow you to gain more insight from your differential expression analysis. First, you will learn to use the biomaRt R package to extract various types of data on your genes of interest. These include translating gene IDs into gene names and vice versa, extracting specific information on genes such as chromosome location, and mapping genes to their homologues in a different organism (for example, mouse to human). You will also learn to subset and extract genes for functional enrichment analysis, both using the R package clusterProfiler as well as other online tools such as Panther, GOrilla.
<br>

### REQUIREMENTS TO PARTICIPATE IN THIS SESSION:
 
Make sure you have the following packages installed in Rstudio: biomaRt, clusterProfiler (note they are both Bioconductor packages).  

Example of how to install Bioconductor R packages:
```r
### Installing biomaRt:
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("biomaRt")

# if the installation worked, you should be able to run the following line without error messages:
library(biomaRt)
```

### OFFICE HOURS:
 
Friday Dec. 4th
* 11am to noon  
* 1:30pm until the start of the session at 2:30pm  
 
Zoom link:
https://nyulangone.zoom.us/j/98456141720
