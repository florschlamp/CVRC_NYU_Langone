## 
## Friday November 6th
## Day 4: “Data Analysis, Part 1: The basics of Differential Expression analysis with DESeq2”
##

## RNA-seq data analysis using DESeq2 and public 'airway' test data
# ===
# RNA-Seq experiment on four human airway smooth muscle cell lines 
# treated with dexamethasone
# ===

## before session: have DESeq2 installed

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("DESeq2")

## Step 1: call ("activate") the library/package we are going to use
library(DESeq2)


## Step 2: set working directory to source file (this script) location
getwd()
# Session > Set Working Directory > To Source File Location


## Step 3: read in data files
# counts and metadata
load("StartData.Rdata")

head(raw_counts_filt)
sample_table

# genes that have 0 counts across all samples already removed (Day 2)
dim(raw_counts_filt)
# rows x columns
table(rowSums(raw_counts_filt) != 0)


## Step 4: load data as DESeq2 object (EXPLORATION) ##

# you should have the same number of raw counts columns and metadata rows
dim(raw_counts_filt)
dim(sample_table)

ddsMat <- DESeqDataSetFromMatrix(countData = raw_counts_filt,
                                 colData = sample_table,
                                 design = ~ 1) # no design needed YET

rld <- rlog(ddsMat) # takes some time to run
# note: if you have too many samples, might want to use vst() instead of rlog()

## basic PCA (integrated DESeq2 plotting function)
# double check metadata for grouping options
sample_table

plotPCA(rld, intgroup=c("dex"))

# or, if you care about the separation between cell type:
plotPCA(rld, intgroup=c("celltype"))

# or both
plotPCA(rld, intgroup=c("dex","celltype"))
# (not very useful if no replicates)

## better PCA plot (extract data, use ggplot)
PCA_data <- plotPCA(rld, intgroup = c("dex","celltype"), returnData=TRUE)
percentVar <- round(100 * attr(PCA_data, "percentVar"))

library(ggplot2)

ggplot(PCA_data, aes(PC1, PC2, color=dex)) + 
  geom_point(aes(shape=celltype), size=4) + 
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))

# maybe flip shape/color ?
ggplot(PCA_data, aes(PC1, PC2, color=celltype)) + 
  geom_point(aes(shape=dex), size=4) + 
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))


# final form:
ggplot(PCA_data, aes(PC1, PC2, color=dex)) + 
  geom_point(aes(shape=celltype), size=4) + 
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme(axis.text=element_text(size=16),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16),
        axis.title.x=element_text(size=18, margin=margin(20,0,0,0)),
        axis.title.y=element_text(size=18, margin=margin(0,20,0,0)))
# maybe add a plot title?


## Step 5: load data as DESeq2 object (ANALYSIS) ##

ddsMat <- DESeqDataSetFromMatrix(countData = raw_counts_filt,
                                 colData = sample_table,
                                 design = ~ dex)
# "design" depends on comparison you want to make!
# what is the biological question you are trying to answer?

ddsMat = DESeq(ddsMat) # normalization and differential expression steps are done here
resultsNames(ddsMat) # print result terms

res = results(ddsMat, name = "dex_treated_vs_control") # extract results for desired term
res # <--- save this table for DE results

# how many DE genes do we have?
# pick significance threshold to assign 'DE'
table(res$padj < 0.05)
# 2175 genes have a significance of adjusted p-value < 0.05

table(res$pvalue < 0.05)
# 3875 genes have a significance of p-value < 0.05

# which are the top DE genes?
head(res[order(res$padj),],10) # global ranking by adj p-value


# explore other cutoffs (logFC)

table(res$padj < 0.05 & res$log2FoldChange > 2)
# 99 genes have significance at padj < 0.05 AND are upregulated by more than 2 logFC

table(res$padj < 0.05 & res$log2FoldChange < -2)
# 70 genes have significance at padj < 0.05 AND are downregulated by move than 2 logFC

table(res$padj < 0.05 & abs(res$log2FoldChange) > 2)
# 169 genes have significance of padj < 0.05 AND are up- or down-regulated by more than 2 logFC
# 169 = 99 + 70 (sanity check)


## save table of results for later
results_dex_treated_vs_control <- data.frame(res)

head(results_dex_treated_vs_control) # save as csv
head(res)

top_10 <- head(results_dex_treated_vs_control[order(results_dex_treated_vs_control$padj),],10)

data_to_plot <- data.frame(genes = rownames(top_10),
                           logFC = top_10$log2FoldChange)

# plot log fold change for each gene
ggplot(data_to_plot, aes(x=genes, y=logFC)) +
  geom_bar(stat="identity")

# flip axis names
ggplot(data_to_plot, aes(x=genes, y=logFC)) +
  geom_bar(stat="identity") +
  theme(axis.text.x = element_text(angle = 90))

# adding error bars
data_to_plot <- data.frame(genes = rownames(top_10),
                           logFC = top_10$log2FoldChange,
                           SE = top_10$lfcSE)

ggplot(data_to_plot, aes(x=genes, y=logFC)) +
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin=logFC-SE, ymax=logFC+SE)) +
  theme(axis.text.x = element_text(angle = 90))

# improved plot:
ggplot(data_to_plot, aes(x=genes, y=logFC)) +
  geom_bar(stat="identity",color="black",fill="#fb8072",width=0.6) +
  geom_errorbar(aes(ymin=logFC-SE, ymax=logFC+SE), width=.2) +
  geom_hline(yintercept=0) +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(title="Top 10 most significant DE genes",
       subtitle="dex treatment vs. control")

## more colors and palettes:
# https://colorbrewer2.org/#type=sequential&scheme=BuGn&n=3 



# extract norm counts
counts_data_norm <- counts(ddsMat, normalized=TRUE)
counts_data_norm_log2 <- log2(1 + counts_data_norm) # +1 is a small pseudocount addition
dim(counts_data_norm_log2)
head(counts_data_norm_log2) # usually between 0 and 20

# save data for later
save(results_dex_treated_vs_control,counts_data_norm_log2,file="DESeq_results.Rdata")

head(results_dex_treated_vs_control)
head(counts_data_norm_log2)

full_results <- cbind(results_dex_treated_vs_control,counts_data_norm_log2)

write.csv(full_results, file="full_results_dex_treated_vs_control.csv")



## for more information on DESeq2 package:
# https://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html


## Example of different designs (~dex + celltype)

ddsMat <- DESeqDataSetFromMatrix(countData = raw_counts_filt,
                                 colData = sample_table,
                                 design = ~ dex + celltype)

ddsMat = DESeq(ddsMat) 
resultsNames(ddsMat) 

res = results(ddsMat, name = "dex_treated_vs_control") # extract results for desired term
res # <--- save this table for DE results

# how many DE genes do we have?
# pick significance threshold to assign 'DE'
table(res$padj < 0.05)
# 3323 genes (instead of 2175 genes) have a significance of adjusted p-value < 0.05

table(res$pvalue < 0.05)
# 4780 genes (instead of 3875 genes) have a significance of p-value < 0.05
