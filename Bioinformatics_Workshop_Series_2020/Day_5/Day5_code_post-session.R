## 
## Friday November 20th
## Day 5: “Data Analysis, Part 2: Publication-Ready RNAseq Plots”
##

## Visualization of RNA-seq data analysis using public 'airway' test data
# ===
# RNA-Seq experiment on four human airway smooth muscle cell lines 
# treated with dexamethasone
# ===

## before session:
# have 'ggplot2', 'ggrepel', 'plotly', 'pheatmap' installed

# code:
# install.packages("ggrepel")
# install.packages("plotly")

## 1) set working directory to source file (this script) location
getwd()
# Session > Set Working Directory > To Source File Location

## 2) read in data files
# counts and metadata
load("StartData.Rdata")
load("DESeq_results.Rdata")

sample_table # metadata table
head(raw_counts_filt) # raw counts
head(counts_data_norm_log2) # normalized counts
head(results_dex_treated_vs_control) # DE results

## Let's start with volcano plots!
library(ggplot2)

# first step: are you going to plot pvalues or adjusted pvalues?
# In general, adjusted pvalues (padj) are better, so I always start there

# second step: get rid of genes that have NAs in the padj column
table(is.na(results_dex_treated_vs_control$pvalue))
table(is.na(results_dex_treated_vs_control$padj))

results <- results_dex_treated_vs_control
results_filt <- results[complete.cases(results$padj),]

dim(results)
dim(results_filt)

head(results_filt)

# prepare data for volcano plot:
data_to_plot <- data.frame(gene = rownames(results_filt), 
                           pvalue = -log10(results_filt$padj), # -log10 scaling 
                           logFC = results_filt$log2FoldChange)
head(data_to_plot) # <- this is what will be plotted

# most basic volcano plot
ggplot(data_to_plot, aes(x = logFC, y = pvalue)) +
  geom_point()

# label it!
ggplot(data_to_plot, aes(x = logFC, y = pvalue)) +
  geom_point() +
  ggtitle(label = "Volcano Plot", subtitle= "dex treatment vs control") +  
  xlab("log2 (Fold Change)") + 
  ylab("-log10 (adjusted p-value)")

# scale the y-axis
ggplot(data_to_plot, aes(x = logFC, y = pvalue)) +
  geom_point() +
  ggtitle(label = "Volcano Plot", subtitle= "dex treatment vs control") +  
  xlab("log2 (Fold Change)") + 
  ylab("-log10 (adjusted p-value)") +
  scale_y_continuous(trans = "log1p")

# general improvements
ggplot(data_to_plot, aes(x = logFC, y = pvalue)) +
  geom_point(size = 1.25, alpha = 0.5) + # alpha is for transparency
  theme_bw(base_size = 16) + 
  ggtitle(label = "Volcano Plot", subtitle= "dex treatment vs control") +  
  xlab("log2 (Fold Change)") + 
  ylab("-log10 (adjusted p-value)") +
  scale_y_continuous(trans = "log1p")

# lets add guiding lines
ggplot(data_to_plot, aes(x = logFC, y = pvalue)) +
  geom_point(size = 1.25, alpha = 0.5) +
  theme_bw(base_size = 16) + 
  ggtitle(label = "Volcano Plot", subtitle= "dex treatment vs control") +  
  xlab("log2 (Fold Change)") + 
  ylab("-log10 (adjusted p-value)") +
  scale_y_continuous(trans = "log1p") +
  geom_vline(xintercept = 0) +   # vertical line around 0
  geom_hline(yintercept = -log10(0.05))  # horizontal line at pvalue cutoff

## more elaborate additions:
# color code genes by expression

pval_cutoff = -log10(0.05)

data_to_plot$color <- ifelse(data_to_plot$logFC > 0 & data_to_plot$pvalue > pval_cutoff,
                             yes = "upregulated",
                             no = ifelse(data_to_plot$logFC < 0 & data_to_plot$pvalue > pval_cutoff,
                                         yes = "downregulated",
                                         no = "not-significant"))

head(data_to_plot)

ggplot(data_to_plot, aes(x = logFC, y = pvalue, color=color)) + # add color variable
  geom_point(size = 1.25, alpha = 0.5) +
  theme_bw(base_size = 16) + 
  ggtitle(label = "Volcano Plot", subtitle= "dex treatment vs control") +  
  xlab("log2 (Fold Change)") + 
  ylab("-log10 (adjusted p-value)") +
  scale_y_continuous(trans = "log1p") +
  geom_vline(xintercept = 0) +   
  geom_hline(yintercept = pval_cutoff) +
  scale_color_manual(values = c("upregulated" = "red", # add manual colors for each group
                                "downregulated" = "blue",
                                "not-significant" = "grey"))

# if you want to remove the legend:
ggplot(data_to_plot, aes(x = logFC, y = pvalue, color=color)) +
  geom_point(size = 1.25, alpha = 0.5) +
  theme_bw(base_size = 16) + ## keep main theme call first
  theme(legend.position = "none") + # remove legend <---
  ggtitle(label = "Volcano Plot", subtitle= "dex treatment vs control") +  
  xlab("log2 (Fold Change)") + 
  ylab("-log10 (adjusted p-value)") +
  scale_y_continuous(trans = "log1p") +
  geom_vline(xintercept = 0) +   
  geom_hline(yintercept = pval_cutoff) +
  scale_color_manual(values = c("upregulated" = "red", 
                                "downregulated" = "blue",
                                "not-significant" = "grey"))


##### my final preferred volcano plot style:
ggplot(data_to_plot, aes(x = logFC, y = pvalue, color=color)) +
  geom_point(size = 1.25, alpha = 0.5) +
  theme_bw(base_size = 16) +
  theme(legend.position = "none") + 
  ggtitle(label = "Volcano Plot", subtitle= "dex treatment vs control") +  
  xlab("log2 (Fold Change)") + 
  ylab("-log10 (adjusted p-value)") +
  scale_y_continuous(trans = "log1p") +
  geom_vline(xintercept = 0) +   
  geom_hline(yintercept = pval_cutoff) +
  scale_color_manual(values = c("upregulated" = "#E64B35",
                                "downregulated" = "#3182bd",
                                "not-significant" = "#636363"))

## a crowd favorite, labeling genes of interest
library(ggrepel)

genes_of_interest <- c("ENSG00000167641","ENSG00000131242","ENSG00000001084") 

subset <- data_to_plot[data_to_plot$gene %in% genes_of_interest,]

head(data_to_plot)
data_to_plot[data_to_plot$gene %in% genes_of_interest,]$color <- "highlight"
data_to_plot[data_to_plot$gene %in% genes_of_interest,]

head(subset)

ggplot(data_to_plot, aes(x = logFC, y = pvalue, color=color)) +
  geom_point(size = 1.25, alpha = 0.5) +
  theme_bw(base_size = 16) +
  theme(legend.position = "none") + 
  ggtitle(label = "Volcano Plot", subtitle= "dex treatment vs control") +  
  xlab("log2 (Fold Change)") + 
  ylab("-log10 (adjusted p-value)") +
  scale_y_continuous(trans = "log1p") +
  geom_vline(xintercept = 0) +   
  geom_hline(yintercept = pval_cutoff) +
  scale_color_manual(values = c("upregulated" = "#E64B35",
                                "downregulated" = "#3182bd",
                                "not-significant" = "#636363",
                                "highlight" = "black")) +
  geom_text_repel(data = subset, mapping = aes(label = gene),
                  size = 5, color = "black",
                  fontface = 'bold',
                  box.padding = unit(0.5, "lines"),
                  point.padding = unit(0.5, "lines"))


### how do we save plots?
ggplot(data_to_plot, aes(x = logFC, y = pvalue, color=color)) +
  geom_point(size = 1.25, alpha = 0.5) +
  theme_bw(base_size = 16) +
  theme(legend.position = "none") + 
  ggtitle(label = "Volcano Plot", subtitle= "dex treatment vs control") +  
  xlab("log2 (Fold Change)") + 
  ylab("-log10 (adjusted p-value)") +
  scale_y_continuous(trans = "log1p") +
  geom_vline(xintercept = 0) +   
  geom_hline(yintercept = pval_cutoff) +
  scale_color_manual(values = c("upregulated" = "#E64B35",
                                "downregulated" = "#3182bd",
                                "not-significant" = "#636363",
                                "highlight" = "black")) +
  geom_text_repel(data = subset, mapping = aes(label = gene),
                  size = 5, color = "black",
                  fontface = 'bold',
                  box.padding = unit(0.5, "lines"),
                  point.padding = unit(0.5, "lines")) +
  ggsave("volcano_plot_genes_labeled_padj0.05.pdf", width=6, height=6)
# you can also save as .png or .jpeg

## another crowd favorite, interactive plot!
library(plotly)

data_to_plot$color <- ifelse(data_to_plot$logFC > 0 & 
                               data_to_plot$pvalue > pval_cutoff,
                             yes = "upregulated",
                             no = ifelse(data_to_plot$logFC < 0 & 
                                           data_to_plot$pvalue > pval_cutoff,
                                         yes = "downregulated",
                                         no = "not-significant"))

head(data_to_plot)

# run plot again but save as variable instead of plotting
p1 <- ggplot(data_to_plot, aes(x = logFC, y = pvalue, color=color, text=gene)) +
  geom_point(size = 1.25, alpha = 0.5) +
  theme_bw(base_size = 16) +
  theme(legend.position = "none") + 
  ggtitle(label = "Volcano Plot", subtitle= "dex treatment vs control") +  
  xlab("log2 (Fold Change)") + 
  ylab("-log10 (adjusted p-value)") +
  scale_y_continuous(trans = "log1p") +
  geom_vline(xintercept = 0) +   
  geom_hline(yintercept = pval_cutoff) +
  scale_color_manual(values = c("upregulated" = "#E64B35",
                                "downregulated" = "#3182bd",
                                "not-significant" = "#636363"))

ggplotly(p1)

# save plot as interactive HTML
htmlwidgets::saveWidget(as_widget(ggplotly(p1)), "volcano_plot_interactive.html")


##
##
#### MA plot
# started being used with microarrays
# visualizes the differences between measurements taken in two samples, 
# by plotting log fold change (M) vs mean expression (average of normalized counts) (A)

data_to_plot <- data.frame(gene = rownames(results_filt), 
                           logFC = results_filt$log2FoldChange,
                           baseMean = results_filt$baseMean)
head(data_to_plot)

ggplot(data_to_plot, aes(x = baseMean, y = logFC)) +
  geom_point(size = 1.25, alpha = 0.5) +
  theme_bw(base_size = 16) +
  theme(legend.position = "none") + 
  ggtitle(label = "MA Plot", subtitle= "dex treatment vs control") +  
  xlab("baseMean") + 
  ylab("log2 (Fold Change)") +
  scale_x_log10()


# color by significance
data_to_plot <- data.frame(gene = rownames(results), 
                           logFC = results$log2FoldChange,
                           baseMean = results$baseMean,
                           pvalue = results$padj)

head(data_to_plot)

data_to_plot$color <- ifelse(data_to_plot$pvalue < 0.05,
                             yes = "sig", no = "non-sig")

head(data_to_plot)

ggplot(data_to_plot, aes(x = baseMean, y = logFC, color=color)) +
  geom_point(size = 1.25, alpha = 0.5) +
  theme_bw(base_size = 16) +
  #theme(legend.position = "none") + 
  ggtitle(label = "MA Plot", subtitle= "dex treatment vs control") +  
  xlab("baseMean") + 
  ylab("log2 (Fold Change)") +
  scale_x_log10()

# genes with low counts don't even get adjusted p-values (NAs, in grey)
# show difference when using p-values

# add logFC thresholds

ggplot(data_to_plot, aes(x = baseMean, y = logFC, color=color)) +
  geom_point(size = 1.25, alpha = 0.5) +
  theme_bw(base_size = 16) +
  #theme(legend.position = "none") + 
  ggtitle(label = "MA Plot", subtitle= "dex treatment vs control") +  
  xlab("baseMean") + 
  ylab("log2 (Fold Change)") +
  scale_x_log10() + 
  geom_hline(yintercept = 2, linetype = 'dashed') + 
  geom_hline(yintercept = -2, linetype = 'dashed')
# try adding a vertical line at baseMean = 30

ggplot(data_to_plot, aes(x = baseMean, y = logFC, color=color)) +
  geom_point(size = 1.25, alpha = 0.5) +
  theme_bw(base_size = 16) +
  #theme(legend.position = "none") + 
  ggtitle(label = "MA Plot", subtitle= "dex treatment vs control") +  
  xlab("baseMean") + 
  ylab("log2 (Fold Change)") +
  scale_x_log10() + 
  geom_hline(yintercept = 2, linetype = 'dashed') + 
  geom_hline(yintercept = -2, linetype = 'dashed') +
  geom_vline(xintercept = 30)

##
##
#### heatmaps
library(pheatmap)

## pick most interesting genes
### cut by pval, cut by logfc
subset_results <- results[results$padj < 0.05 &
                         abs(results$log2FoldChange) > 2 &
                         results$baseMean > 30,]
dim(subset_results)
# 133 genes

subset_norm_counts <- counts_data_norm_log2[rownames(counts_data_norm_log2) %in%
                                              rownames(subset_results),]
dim(subset_norm_counts)

# most simple heatmap
pheatmap(subset_norm_counts)

pheatmap(subset_norm_counts,
         show_rownames = F)

pheatmap(subset_norm_counts,
         show_rownames = F,
         cellwidth = 20,
         angle_col = 0)

head(subset_norm_counts)

# scale by row
pheatmap(subset_norm_counts,
         show_rownames = F,
         cellwidth = 20,
         angle_col = 0,
         scale = "row")


## WARNING ##
# rank-Z scaled
# interpretation issues #

test_data <- data.frame(row.names = c("gene1","gene2","gene3"),
                        S1 = c(1,10,0.1),
                        S2 = c(2,20,0.3),
                        S3 = c(3,30,0.5))
# not scaled
pheatmap(test_data,
         cluster_rows = F,
         cluster_cols = F,
         cellwidth = 40,
         cellheight = 30,
         display_numbers = test_data,
         angle_col = 0)

# scaled
pheatmap(test_data,
         cluster_rows = F,
         cluster_cols = F,
         cellwidth = 40,
         cellheight = 30,
         display_numbers = test_data,
         angle_col = 0,
         scale = "row")

## END OF WARNING ##

# annotation key rownames must be same as data colnames


annotation_key <- data.frame(row.names = sample_table$sample_name,
                             treatment = sample_table$dex)

colnames(subset_norm_counts)
rownames(annotation_key)

pheatmap(subset_norm_counts,
         show_rownames = F,
         cellwidth = 20,
         angle_col = 0,
         annotation_col = annotation_key)

# what if we care about cell type too?
# (not in this case, since we did the DE analysis for treatment alone)
annotation_key <- data.frame(row.names = sample_table$sample_name,
                             treatment = sample_table$dex,
                             celltype = sample_table$celltype)

pheatmap(subset_norm_counts,
         show_rownames = F,
         cellwidth = 20,
         angle_col = 0,
         annotation_col = annotation_key)

