## 
## Friday October 23rd
## Day 3: “Data Exploration, Part 2: Basic RNA-seq Plots”
##

# before session:
# have 'ggplot2', 'pheatmap', 'RColorBrewer', 'reshape2' installed

# code:
# install.packages("ggplot2")


## Step 1: call ("activate") libraries (packages) we are going to use
library(ggplot2)

## Step 2: set working directory to source file (this script) location
getwd()
# Session > Set Working Directory > To Source File Location


## Step 3: read in data files
# counts and metadata
load("StartData.Rdata")

head(raw_counts_filt)
sample_table

## Step 4: normalize counts by library size and log scale
lib_sizes=colSums(raw_counts_filt)
size.factor=lib_sizes/exp(mean(log(lib_sizes)))
norm_counts=t(t(raw_counts_filt)/size.factor)
head(norm_counts)

norm_log_counts = log2(norm_counts+1)
head(norm_log_counts)

# Plot 1) comparing replicates
# pairwise scatterplot

norm_log_counts <- data.frame(norm_log_counts)

## base plot
# basic
plot(norm_log_counts$C1,norm_log_counts$C2)

## ggplot
# basic
ggplot(norm_log_counts, aes(x=C1, y=C2)) +
  geom_point() 

# add line
ggplot(norm_log_counts, aes(x=C1, y=C2)) + 
  geom_point() + geom_abline(slope=1, color="red")

# add better labels
ggplot(norm_log_counts, aes(x=C1, y=C2)) +
  geom_point() + geom_abline(slope=1, color="red") +
  labs(title="pairwise scatterplot between samples C1 and C2",
       x="norm counts for sample C1",
       y="norm counts for sample C2")

# themes
ggplot(norm_log_counts, aes(x=C1, y=C2)) +
  geom_point() + theme_bw()

ggplot(norm_log_counts, aes(x=C1, y=C2)) +
  geom_point() + theme_classic()

ggplot(norm_log_counts, aes(x=C1, y=C2)) +
  geom_point() + theme_grey() # default?

ggplot(norm_log_counts, aes(x=C1, y=C2)) +
  geom_point() + theme_minimal()
  
  
# manipulate points
ggplot(norm_log_counts, aes(x=C1, y=C2)) +
  geom_point(alpha=0.25) 

ggplot(norm_log_counts, aes(x=C1, y=C2)) +
  geom_point(alpha=0.25, color="blue") 




# filtered counts
norm_log_counts_filt <- read.csv("normalized_counts_after_log_and_filtering.csv", row.names = 1)

ggplot(norm_log_counts_filt, aes(x=C1, y=C2)) +
  geom_point(alpha=0.25, color="blue") 

# let's add some stats
rsq <- cor(norm_log_counts_filt$C1,norm_log_counts_filt$C2) ^ 2
rsq

ggplot(norm_log_counts_filt, aes(x=C1, y=C2)) +
  geom_point(alpha=0.25) + geom_abline(slope=1, color="red") +
  labs(title="pairwise scatterplot between samples C1 and C2",
       subtitle="R squared = 0.8649412",
       x="norm counts for sample C1",
       y="norm counts for sample C2")

ggplot(norm_log_counts_filt, aes(x=C1, y=C2)) +
  geom_point(alpha=0.25) + geom_abline(slope=1, color="red") +
  labs(title="pairwise scatterplot between samples C1 and C2",
       subtitle=paste0("R squared = ",rsq),
       x="norm counts for sample C1",
       y="norm counts for sample C2")

round(rsq,digits=3)

# final form:
ggplot(norm_log_counts_filt, aes(x=C1, y=C2)) +
  geom_point(alpha=0.25) + geom_abline(slope=1, color="red") +
  labs(title="pairwise scatterplot between samples C1 and C2",
       subtitle=paste0("R squared = ",round(rsq,digits=3)),
       x="norm counts for sample C1",
       y="norm counts for sample C2")



## Any questions so far?



cor(norm_log_counts_filt$C1,norm_log_counts_filt$C2) ^ 2
cor(norm_log_counts_filt$C1,norm_log_counts_filt$T1) ^ 2


# Plot 2) Comparing multiple samples
# heatmap of sample distances
library(pheatmap)

# calculate distances
sampleDists <- dist(t(norm_log_counts_filt))
sampleDistMatrix <- as.matrix(sampleDists)

# simple heatmap
pheatmap(sampleDistMatrix)

# change colors manually
library(RColorBrewer)
## for more colors: https://colorbrewer2.org/#type=sequential&scheme=BuGn&n=3

colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(50)

pheatmap(sampleDistMatrix,
         col=colors)

# add group labels
sample_table
rownames(sampleDistMatrix)

annotation <- data.frame(row.names = sample_table$sample_name, #HAS to be the same as matrix
                         cell_type = sample_table$celltype,
                         treatment = sample_table$dex)

annotation

pheatmap(sampleDistMatrix,
         col=colors,
         annotation_col = annotation)



## show examples from other datasets



# Plot 3) Plotting genes of interest
# A) plot individual genes

gene_of_interest <- "ENSG00000179094"

norm_counts_subset <- norm_log_counts[gene_of_interest,]
norm_counts_subset

data_to_plot <- data.frame(row.names = sample_table$sample_name,
                           cell_type = sample_table$celltype,
                           treatment = sample_table$dex,
                           norm_counts = as.numeric(norm_counts_subset))

data_to_plot

# basic plot
ggplot(data_to_plot, aes(x=treatment, y=norm_counts)) +
  geom_point()

ggplot(data_to_plot, aes(x=cell_type, y=norm_counts)) +
  geom_point()

ggplot(data_to_plot, aes(x=treatment, y=norm_counts)) +
  geom_point(size=3,aes(color=cell_type))

# add lines
ggplot(data_to_plot, aes(x=treatment, y=norm_counts, group=cell_type)) +
  geom_point(size=3,aes(color=cell_type)) +
  geom_line()

ggplot(data_to_plot, aes(x=treatment, y=norm_counts, group=cell_type)) +
  geom_point(size=3,aes(color=cell_type)) +
  geom_line(aes(color=cell_type))

# boxplots instead
ggplot(data_to_plot, aes(x=treatment, y=norm_counts)) +
  geom_boxplot()

ggplot(data_to_plot, aes(x=treatment, y=norm_counts)) +
  geom_boxplot(aes(fill=treatment))

# add points
ggplot(data_to_plot, aes(x=treatment, y=norm_counts)) +
  geom_boxplot() + geom_point()

# add labels
ggplot(data_to_plot, aes(x=treatment, y=norm_counts)) +
  geom_boxplot(aes(fill=treatment)) + geom_point() +
  labs(title=gene_of_interest)


# B) plot multiple genes
genes_of_interest <- c("ENSG00000152583",
                       "ENSG00000179094",
                       "ENSG00000116584",
                       "ENSG00000189221",
                       "ENSG00000120129",
                       "ENSG00000148175",
                       "ENSG00000178695",
                       "ENSG00000109906",
                       "ENSG00000134686",
                       "ENSG00000101347")

norm_counts_subset <- norm_log_counts[genes_of_interest,]
norm_counts_subset

# flip the table
t(norm_counts_subset)
sample_table

metadata_subset <- data.frame(row.names = sample_table$sample_name,
                              cell_type = sample_table$celltype,
                              treatment = sample_table$dex)
metadata_subset

data_to_plot <- cbind(metadata_subset,t(norm_counts_subset))
data_to_plot

# 'melt' the table
library(reshape2)

data_to_plot_melted <- melt(data_to_plot)
head(data_to_plot_melted)

colnames(data_to_plot_melted) <- c("cell_type","treatment","gene","norm_counts")
head(data_to_plot_melted)

# plot by gene
ggplot(data_to_plot_melted, aes(x=gene, y=norm_counts)) +
  geom_point()

ggplot(data_to_plot_melted, aes(x=gene, y=norm_counts)) +
  geom_point() +
  theme(axis.text.x = element_text(angle = 90))

ggplot(data_to_plot_melted, aes(x=gene, y=norm_counts)) +
  geom_point(aes(color=treatment)) +
  theme(axis.text.x = element_text(angle = 90))

# plot by treatment
ggplot(data_to_plot_melted, aes(x=treatment, y=norm_counts)) +
  geom_point(aes(color=gene))

ggplot(data_to_plot_melted, aes(x=treatment, y=norm_counts)) +
  geom_point(aes(color=gene)) +
  stat_summary(fun = mean, geom="line", mapping = aes(group=gene, color=gene))

# final form:
ggplot(data_to_plot_melted, aes(x=treatment, y=norm_counts)) +
  geom_point(aes(color=gene)) +
  stat_summary(fun = mean, geom="line", mapping = aes(group=gene, color=gene)) +
  labs(title="normalized counts of top 10 most variable genes",
       subtitle="grouped by treatment",
       y="log transformed normalized counts",
       x="dex status")


# adding more stats:
# smaller subset (10 genes is too much)

genes_of_interest <- c("ENSG00000152583",
                       "ENSG00000179094",
                       "ENSG00000116584")

norm_counts_subset <- norm_log_counts[genes_of_interest,]
rownames(norm_counts_subset) <- c("SPARCL1","PER1","ARHGEF2")
metadata_subset <- data.frame(row.names = sample_table$sample_name,
                              cell_type = sample_table$celltype,
                              treatment = sample_table$dex)
data_to_plot <- cbind(metadata_subset,t(norm_counts_subset))
data_to_plot_melted <- melt(data_to_plot)
colnames(data_to_plot_melted) <- c("cell_type","treatment","gene","norm_counts")
head(data_to_plot_melted)

# plot by gene
ggplot(data_to_plot_melted, aes(x=gene, y=norm_counts)) +
  geom_point(aes(color=treatment))

# boxplot?
ggplot(data_to_plot_melted, aes(x=gene, y=norm_counts)) +
  geom_boxplot(aes(fill=treatment))

# points might be better
ggplot(data_to_plot_melted, aes(x=gene, y=norm_counts)) +
  geom_point(aes(color=treatment)) +
  stat_summary(fun.min=function(x)(mean(x)-sd(x)), 
               fun.max=function(x)(mean(x)+sd(x)),
               geom="errorbar", aes(group=treatment),
               width=0.1, position=position_dodge(.4)) +
  stat_summary(fun=mean, geom="point", shape=23, color="black", aes(fill=treatment), 
               size=2,position=position_dodge(.4)) 



## Homework for next session:
# install the DESeq2 package (from Bioconductor)




