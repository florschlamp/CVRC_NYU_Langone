## 
## Friday October 9th
## Day 2: “Data Exploration, Part 1: Basic RNA-seq Data Manipulation”
##

## Manual RNA-seq data exploration using public 'airway' test data
# ===
# RNA-Seq experiment on four human airway smooth muscle cell lines 
# treated with dexamethasone
# ===

# before session:
# have ggplot2 and edgeR installed

# code:
# install.packages("ggplot2")
# 
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("edgeR")



## Step 1: call ("activate") libraries (packages) we are going to use
library(ggplot2)
library(edgeR)



## Step 2: set working directory to source file (this script) location
getwd()
# Session > Set Working Directory > To Source File Location



## Step 3: read in data files
# counts:
raw_counts <- read.csv("airway_scaledcounts.csv",row.names = 1)
head(raw_counts)
dim(raw_counts) # 38,694 genes (rows) by 8 samples (columns)

# metadata:
sample_table <- read.csv("airway_metadata.csv")
sample_table
dim(sample_table) # 8 samples (rows) with 4 info columns

# number of columns (samples) in raw counts should be equal to number of rows in metadata


## Step 4: metadata manipulation
# we don't care about the geo_id right now
dim(sample_table)
sample_table$geo_id <- NULL

dim(sample_table)

sample_table




# lets use more informative sample names
better_sample_names <- c("C1","T1","C2","T2","C3","T3","C4","T4")
sample_table$sample_name <- better_sample_names
sample_table

head(raw_counts)
colnames(raw_counts)
sample_table$id
colnames(raw_counts) <- sample_table$sample_name

head(raw_counts)




## Step 5: explore library size
lib_sizes = colSums(raw_counts)
lib_sizes
summary(lib_sizes)
summary(lib_sizes/1000000) # in million reads
summary(lib_sizes)/1000000

data_to_plot <- data.frame(sample_name = sample_table$sample_name,
                           libsize = lib_sizes)
data_to_plot

data_to_plot$sample_name


# basic ggplot (bar plot)
ggplot(data_to_plot, aes(x=sample_name, y=lib_sizes)) +
  geom_bar(stat="identity")

# improve plot a bit
ggplot(data_to_plot, aes(x=sample_name, y=lib_sizes/1000000)) +
  geom_bar(stat="identity") + 
  labs(y="Library size (in million of reads)", 
       title="Library size per sample (in million of reads)")



# improve more
summary(lib_sizes/1000000) # mean is 22.6

ggplot(data_to_plot, aes(x = reorder(sample_name, -lib_sizes/1000000), y=lib_sizes/1000000)) +
  geom_bar(stat="identity") + 
  coord_flip() +
  labs(y="Library size (in million of reads)", 
       title="Library size per sample (in million of reads)") +
  geom_hline(yintercept=22.6)

# 'hard code' of mean value is not ideal
# different ways of calculating the mean:
mean(lib_sizes/1000000)
mean(lib_sizes)/1000000
summary(lib_sizes/1000000)[[4]]
as.numeric(summary(lib_sizes/1000000)["Mean"])


ggplot(data_to_plot, aes(x = reorder(sample_name, -lib_sizes/1000000), y=lib_sizes/1000000)) +
  geom_bar(stat="identity") + 
  coord_flip() +
  labs(y="Library size (in million of reads)", title="Library size per sample (in million of reads)") +
  geom_hline(yintercept=mean(lib_sizes/1000000)) 


# outliers could be identified in this step already



## Step 6: remove genes that have 0 counts across all samples
head(raw_counts)


dim(raw_counts)
table(rowSums(raw_counts) != 0)
genes_to_keep = rowSums(raw_counts) != 0
raw_counts_filt=raw_counts[genes_to_keep,] # when subsetting, TABLE[rows,columns]

dim(raw_counts_filt)
# we went down from 38,694 genes to 25,258 genes

colSums(raw_counts)/1000000
colSums(raw_counts_filt)/1000000




## Optional: save data for later:
save(raw_counts_filt,sample_table,file="StartData.Rdata")


## BREAK ##


## Step 7: normalize counts by library size
head(raw_counts_filt)

lib_sizes=colSums(raw_counts_filt)
size.factor=lib_sizes/exp(mean(log(lib_sizes)))
norm_counts=t(t(raw_counts_filt)/size.factor)

head(norm_counts)

colSums(norm_counts) # <-- they are now all the same

## Step 8: log scale normalized counts
range(raw_counts_filt)
range(norm_counts)

norm_log_counts = log2(norm_counts+1)

head(norm_log_counts)
range(norm_log_counts)

dim(raw_counts)
dim(raw_counts_filt)
dim(norm_counts)
dim(norm_log_counts)

## Step 9: visualize data and filter out genes with low read counts

# plot MultiDimensional Scaling plot (similar to PCA)
# plots distances between gene expression profiles
plotMDS(norm_log_counts)

# plot distribution of counts with a histogram
hist(norm_log_counts)

# make slightly better histogram
hist(norm_counts,breaks=30, main="histogram of normalized counts")
hist(norm_log_counts,breaks=30, main="histogram of log normalized counts")

# choose cutoff to filter out low read counts
keep=rowSums((norm_log_counts)>6) >= 2
# rowSums((norm_log_counts)>6) <-- (number of log norm counts) >= 2 <-- (number of samples)
norm_log_counts_filt <- norm_log_counts[keep,]

dim(norm_log_counts)
dim(norm_log_counts_filt)

hist(norm_log_counts_filt,breaks=30, main="histogram of filtered log normalized counts")

plotMDS(norm_log_counts_filt)

## show examples of when this filtering really makes a difference <---

# Visualization bonus:
head(norm_log_counts_filt)
plotMDS(norm_log_counts_filt, col=c("blue","red","blue","red","blue","red","blue","red") )
plotMDS(norm_log_counts_filt, col=c("blue","blue","red","red","green","green","orange","orange") )



## Optional step: removing outliers or subsetting data

colnames(norm_log_counts_filt)
subset_samples <- c("C3","T3","C4","T4")
subset_tissue3and4 <- norm_log_counts_filt[,subset_samples]

dim(norm_log_counts_filt)
dim(subset_tissue3and4) # 12 thousand genes in 4 samples

head(subset_tissue3and4)

plotMDS(subset_tissue3and4)


## Extra: saving processed data as text files vs. Rdata files

write.csv(norm_log_counts_filt,file="normalized_counts_after_log_and_filtering.csv")

save(norm_log_counts_filt,file="normalized_counts_after_log_and_filtering.Rdata")
# also save multiple ones in one Rdata file
save(raw_counts,raw_counts_filt,norm_counts,norm_log_counts,norm_log_counts_filt,
     file="Day2_all_intermediate_tables.Rdata")

load("Day2_all_intermediate_tables.Rdata")


## Homework for next session:
# install the following packages: 'pheatmap', 'RColorBrewer', 'reshape2'


## ===
## one of the sources for more info:
# https://www.bioconductor.org/packages/release/bioc/vignettes/gage/inst/doc/RNA-seqWorkflow.pdf
